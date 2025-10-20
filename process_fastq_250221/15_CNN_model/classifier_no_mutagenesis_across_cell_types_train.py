#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import train_test_split
from torch.utils.data import Subset
import torch
import random
import torch
import torch.nn as nn
import torch.nn.functional as F
import os

# Set random seed for reproducibility
seed = 0
random.seed(seed)
np.random.seed(seed)    
torch.manual_seed(seed)
torch.cuda.manual_seed(seed)
torch.cuda.manual_seed_all(seed)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False

output_dir = "~/lab_projects/splicing/for_kai/no_mutagenesis_data/"
# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)
# Change to output directory
os.chdir(output_dir)

GX = pd.read_csv('~/lab_projects/splicing/for_kai/V4_data_with_offset_major_event_only/241126_gene_expression_all.csv')
PSI = pd.read_csv('~/lab_projects/splicing/for_kai/V4_data_with_offset_major_event_only/250407_measurement_ratio.csv')
Barcode = pd.read_csv('~/lab_projects/splicing/for_kai/V4_data_with_offset_major_event_only/250407_twist_barcodes_with_splice_site_annotation.csv')


GX_dict = GX.groupby("condition").apply(
    lambda x: np.array(x.drop(columns=["Unnamed: 0", "condition"]).values.flatten().tolist())
).to_dict()


# make index as the barcode
Barcode = Barcode.set_index("index_offset")
Barcode = Barcode.drop(columns=["Unnamed: 0"])

# make index as the barcode
PSI = PSI.set_index("index_offset")
PSI = PSI.drop(columns=["Unnamed: 0"])

num_na = PSI.isna().sum().sum()
num_non_na = PSI.size - num_na

print(f"NA: {num_na}")
print(f"~NA: {num_non_na}")


def onehot_encode_sequences(sequences):
    """
    One-hot encodes a list of DNA sequences.

    Args:
        sequences (list of str): A list of DNA sequences, each containing characters 'A', 'G', 'C', 'T'.

    Returns:
        np.ndarray: A 3D numpy array where each sequence is one-hot encoded along the last axis.
    """
    # Define a mapping for bases
    base_mapping = {'A': 0, 'G': 1, 'C': 2, 'T': 3, 'N': 4}
    n_sequences = len(sequences)
    sequence_length = len(sequences[0]) if sequences else 0
    
    # Initialize the one-hot encoded array
    onehot_array = np.zeros((n_sequences, sequence_length, 4), dtype=int)
    
    for i, seq in enumerate(sequences):
        for j, base in enumerate(seq):
            if base in base_mapping:  # Ensure valid base
                onehot_array[i, j, base_mapping[base]] = 1
                
    return onehot_array

Barcode_onehot = onehot_encode_sequences(list(Barcode['padded_sequence'].values))
# Barcode_onehot = np.concatenate((Barcode_onehot, twist_array), axis=-1)
print(Barcode_onehot.shape)
Barcode_dict = dict(zip(Barcode.index, Barcode_onehot))

class PSIDataset(Dataset):
    def __init__(self, PSI_df, Barcode_dict, GX_dict):
        self.samples = []
    
        for barcode in PSI_df.index:
            if barcode not in Barcode_dict:
                continue  
            
            barcode_array = Barcode_dict[barcode] 
            
            for celltype in PSI_df.columns:
                if celltype not in GX_dict:
                    continue 
                
                gx_array = GX_dict[celltype] 
                psi_value = PSI_df.loc[barcode, celltype]
                
                if not np.isnan(psi_value):
                    self.samples.append((celltype, barcode, barcode_array, gx_array, psi_value))
    
    def __len__(self):
        return len(self.samples)
    
    def __getitem__(self, idx):
        celltype, barcode, barcode_array, gx_array, psi_value = self.samples[idx]
        return celltype, barcode, barcode_array, gx_array, psi_value
    
def split_dataset_by_barcode(dataset, test_size=0.2, random_state=42):
    """
    Split the dataset into train and test subsets based on barcode.
    Params:
    - dataset: PSIDataset
    - test_size: ratio
    - random_state
    
    Return:
    - train_dataset
    - test_dataset
    """
    # Extract all barcodes
    barcodes = np.array([sample[1] for sample in dataset.samples])
    
    # Get unique barcodes and split
    unique_barcodes = np.unique(barcodes)
    
    train_barcodes, test_barcodes = train_test_split(
        unique_barcodes, test_size=test_size, random_state=random_state
    )
    
    print(f"Train barcodes: {train_barcodes[:5]}")
    print(f"Test barcodes: {test_barcodes[:5]}") 
    
    train_barcodes_set = set(train_barcodes)
    test_barcodes_set = set(test_barcodes)
    
    # Get indices for train and test
    train_indices = np.where(np.isin(barcodes, list(train_barcodes_set)))[0]
    test_indices = np.where(np.isin(barcodes, list(test_barcodes_set)))[0]
    print(f"Train indices: {len(train_indices)}, Test indices: {len(test_indices)}")
    
    # Construct Subsets
    train_dataset = Subset(dataset, train_indices)
    test_dataset = Subset(dataset, test_indices)
    
    return train_dataset, test_dataset

def split_dataset_by_celltype(dataset, test_size=0.2, random_state=42):
    """
    Split the dataset into train and test subsets based on celltype.
    Params:
    - dataset: PSIDataset
    - test_size: ratio
    - random_state
    
    Return:
    - train_dataset
    - test_dataset
    """
    celltypes = np.array([sample[0] for sample in dataset.samples])
    
    unique_celltypes = np.unique(celltypes)
    train_celltypes, test_celltypes = train_test_split(
        unique_celltypes, test_size=test_size, random_state=random_state
    )
    
    print(f"Train celltypes: {len(train_celltypes)}, Test celltypes: {len(test_celltypes)}")
    
    train_indices = [i for i, celltype in enumerate(celltypes) if celltype in train_celltypes]
    test_indices = [i for i, celltype in enumerate(celltypes) if celltype in test_celltypes]
    
    train_dataset = Subset(dataset, train_indices)
    test_dataset = Subset(dataset, test_indices)
    
    assert len(train_dataset) == len(train_indices), "Train dataset size mismatch!"
    assert len(test_dataset) == len(test_indices), "Test dataset size mismatch!"
    
    return train_dataset, test_dataset

def split_dataset_random(dataset, test_size=0.2, random_state=42):

    dataset_size = len(dataset)
    
    indices = list(range(dataset_size))
    
    train_indices, test_indices = train_test_split(
        indices, test_size=test_size, random_state=random_state
    )
    
    train_dataset = Subset(dataset, train_indices)
    test_dataset = Subset(dataset, test_indices)
    
    return train_dataset, test_dataset



# ## Convolutional Neural Network

class GexFullyConnected(nn.Module):
    def __init__(self):
        super().__init__()

        self.fc1 = nn.Linear(in_features=19221, out_features=1024)
        self.fc2 = nn.Linear(in_features=1024, out_features=256)
        self.fc3 = nn.Linear(in_features=256, out_features=128)
        self.BN1 = nn.BatchNorm1d(1024)
        self.BN2 = nn.BatchNorm1d(256) 
        self.dropout = nn.Dropout(0.2)     

    def forward(self, x):
        x = F.relu(self.BN1(self.fc1(x)))
        x = F.relu(self.BN2(self.fc2(x)))
        x = F.relu(self.fc3(x))
        return x

class ResNetCNN(nn.Module):
    def __init__(self, channels, dropout=0.2):
        super().__init__() 
        
        self.conv1_5 =  nn.Conv1d(channels, channels, kernel_size=5, padding=2)  
        self.conv1_11 = nn.Conv1d(channels, channels, kernel_size=11, padding=5)
        self.conv1_21 = nn.Conv1d(channels, channels, kernel_size=21, padding=10)

        self.batchnorm1 = nn.BatchNorm1d(channels * 3) 
        self.dropout = nn.Dropout(dropout)
        self.conv_merge = nn.Conv1d(channels * 3, channels, kernel_size=1)

        self.conv2 = nn.Conv1d(channels, channels, kernel_size=5, padding=2)
        self.batchnorm2 = nn.BatchNorm1d(channels)

    def forward(self, x):
        residual = x 
        x1 = self.conv1_5(x)
        x2 = self.conv1_11(x)
        x3 = self.conv1_21(x)

        x = torch.cat([x1, x2, x3], dim=1) 

        x = self.dropout(F.relu(self.batchnorm1(x)))
        x = self.conv_merge(x) 
        x = self.dropout(F.relu(self.batchnorm2(self.conv2(x))))
        
        x += residual  
        return x

class SequenceCNN(nn.Module):
    def __init__(self):
        super().__init__()
        
        self.dropout = nn.Dropout(0.2)
        
        self.conv1 = nn.Conv1d(in_channels=4, out_channels=128, kernel_size=1, padding='same')
        self.batchnorm1 = nn.BatchNorm1d(128)
    
        self.resnetCNN1 = nn.Sequential(
            ResNetCNN(128),
            ResNetCNN(128),
            ResNetCNN(128),
        )
        
        self.maxpool = nn.MaxPool1d(kernel_size=2)
        
        self.conv2 = nn.Conv1d(in_channels=128, out_channels=64, kernel_size=1, padding='same')
        self.batchnorm2 = nn.BatchNorm1d(64)
        
        self.resnetCNN2 = nn.Sequential(
            ResNetCNN(64),
            ResNetCNN(64),
            ResNetCNN(64),
        )
        
        self.conv3 = nn.Conv1d(in_channels=64, out_channels=16, kernel_size=1, padding='same')
        self.batchnorm3 = nn.BatchNorm1d(16)
        
        self.resnetCNN3 = nn.Sequential(
            ResNetCNN(16),
            ResNetCNN(16),
            ResNetCNN(16),
        )

        self.conv4 = nn.Conv1d(in_channels=128, out_channels=16, kernel_size=1, padding='same')
        self.maxpool8 = nn.MaxPool1d(kernel_size=8)
        
        # Linear layers
        self.fc1 = nn.Linear(in_features=496, out_features=128)
        
    def forward(self, x):
        x = self.dropout(F.relu(self.batchnorm1(self.conv1(x))))
        residual = x
        x_residual = x
        x = self.resnetCNN1(x)
        x += residual
        x = self.maxpool(x)
        
        x = self.dropout(F.relu(self.batchnorm2(self.conv2(x))))
        residual = x
        x = self.resnetCNN2(x)
        x += residual
        x = self.maxpool(x)
        
        x = self.dropout(F.relu(self.batchnorm3(self.conv3(x))))
        residual = x
        x = self.resnetCNN3(x)
        x += residual
        x = self.maxpool(x)
        
        x_residual = self.conv4(x_residual)
        x_residual = self.maxpool8(x_residual)
        
        x = x + x_residual
        
        x = x.view(x.size(0), -1)
        x = self.fc1(x)
        
        return x

class Predictor_gene_barcode(nn.Module):
    def __init__(self):
        super().__init__()
        
        self.Barcode_model = SequenceCNN()
        self.GX_model = GexFullyConnected()

        self.fc1 = nn.Linear(in_features=128, out_features=64)
        self.bn1 = nn.BatchNorm1d(64)
        self.fc2 = nn.Linear(in_features=64, out_features=1)
        
        self.dropout = nn.Dropout(0.2)
    
    def forward(self, barcode, gx, condition):
        
        z_barcode = self.Barcode_model(barcode)
        
        z_gx = self.GX_model(gx)
        x = z_barcode + z_gx + condition
        x = self.dropout(F.leaky_relu(self.bn1(self.fc1(x))))
        x = self.fc2(x)
        
        return x
    
class Predictor_barcode(nn.Module):
    def __init__(self):
        super().__init__()
        
        self.Barcode_model = SequenceCNN()

        self.fc1 = nn.Linear(in_features=128, out_features=64)
        self.bn1 = nn.BatchNorm1d(64)
        self.fc2 = nn.Linear(in_features=64, out_features=1)
        
        self.dropout = nn.Dropout(0.2)
    
    def forward(self, barcode):
        
        x = self.Barcode_model(barcode)
        x = self.dropout(F.leaky_relu(self.bn1(self.fc1(x))))
        x = self.fc2(x)
        
        return x


# ## Pretrain classifier
def pretrain_barcode_model(device, train_loader, test_loader, epochs, lr, save_path="DC_TMP_model_barcode_PSI_5epoch.pth"):
    model_barcode = Predictor_barcode()
    model_barcode = model_barcode.to(device)

    MSELoss = nn.MSELoss()
    optimizer = torch.optim.Adam(model_barcode.parameters(), lr=lr)

    for epoch in range(epochs):
        
        for i, (_, _, barcode, gx, psi) in enumerate(train_loader):
            
            barcode = barcode.to(device).float()
            barcode = barcode.permute(0, 2, 1)
            psi = psi.to(device).float()
            
            optimizer.zero_grad()
            
            psi_pred = model_barcode(barcode)
            psi_loss = MSELoss(psi_pred.squeeze(), psi)
            
            psi_loss.backward()
            optimizer.step()
            
            if i % 100 == 0:
                print(f"Epoch {epoch}, Batch {i}, psi_loss: {psi_loss.item()}")
                
        model_barcode.eval()
        predict_psi = []
        real_psi = []
        with torch.no_grad():
            psi_loss = 0
            for i, (_, _, barcode, gx, psi) in enumerate(test_loader):
                barcode = barcode.to(device).float()
                barcode = barcode.permute(0, 2, 1)
                psi = psi.to(device).float()
            
                psi_pred = model_barcode(barcode)
                psi_loss += MSELoss(psi_pred.squeeze(), psi).item()
                
                predict_psi.extend(psi_pred.squeeze().cpu().numpy())
                real_psi.extend(psi.cpu().numpy())
                
            psi_loss /= len(test_loader)
            print(f"Epoch {epoch}, Test PSI loss: {psi_loss}")
        
        model_barcode.train()

    torch.save(model_barcode.state_dict(), save_path)
    return model_barcode

# ## Finetune with gene expression

def finetune_with_gene_expression(device, train_loader, test_loader, epochs, lr, pretrained_model_path, save_path="DC_TMP_model_barcode_gene_PSI_5epoch.pth"):
    model_gene_barcode = Predictor_gene_barcode()
    model_gene_barcode = model_gene_barcode.to(device)

    model_barcode = Predictor_barcode()
    model_barcode = model_barcode.to(device)
    model_barcode.load_state_dict(torch.load(pretrained_model_path))
    model_barcode.eval()

    MSELoss = nn.MSELoss()
    optimizer = torch.optim.Adam(model_gene_barcode.parameters(), lr=lr)

    for epoch in range(epochs):
        
        for i, (_, _, barcode, gx, psi) in enumerate(train_loader):
            
            barcode = barcode.to(device).float()
            barcode = barcode.permute(0, 2, 1)
            gx = gx.to(device).float()
            psi = psi.to(device).float()
            
            optimizer.zero_grad()
            
            condition = model_barcode.Barcode_model(barcode).detach()
            psi_pretrain = model_barcode(barcode).detach()
            
            psi_residual = model_gene_barcode(barcode, gx, condition)
            
            psi_loss = MSELoss(psi_residual.squeeze(), psi-psi_pretrain.squeeze())
            
            psi_loss.backward()
            optimizer.step()
            
            if i % 100 == 0:
                print(f"Epoch {epoch}, Batch {i}, psi_loss: {psi_loss.item()}")
                
        model_gene_barcode.eval()
        predict_psi = []
        real_psi = []
        with torch.no_grad():
            psi_loss = 0
            recon_loss = 0
            for i, (_, _, barcode, gx, psi) in enumerate(test_loader):
                barcode = barcode.to(device).float()
                barcode = barcode.permute(0, 2, 1)
                gx = gx.to(device).float()
                psi = psi.to(device).float()
            
                condition = model_barcode.Barcode_model(barcode)
                psi_pretrain = model_barcode(barcode)
            
                psi_residual = model_gene_barcode(barcode, gx, condition)
                psi_loss += MSELoss(psi_residual.squeeze(), psi-psi_pretrain.squeeze()).item()
                
                predict_psi.extend(psi_residual.squeeze().cpu().numpy()+psi_pretrain.squeeze().cpu().numpy())
                real_psi.extend(psi.cpu().numpy())
                
            psi_loss /= len(test_loader)
            print(f"Epoch {epoch}, Test PSI loss: {psi_loss}")
        
        model_gene_barcode.train()
        
    torch.save(model_gene_barcode.state_dict(), save_path)
    return model_gene_barcode


device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
lr = 1e-4
epochs = 5
batch_size = 1024

dataset = PSIDataset(PSI, Barcode_dict, GX_dict)
dataset_size = len(dataset)
print(f"Total samples: {dataset_size}")

# Call the function for each seed
# Set random seed for reproducibility
# seed_set = [0,1,2,3,4]
seed_set = range(10)
for seed in seed_set:
    random.seed(seed)
    np.random.seed(seed)    
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)

    # Create train/test split using current seed
    # train_dataset, test_dataset = split_dataset_random(dataset, test_size=0.2, random_state=seed)
    # train_dataset, test_dataset = split_dataset_by_barcode(dataset, test_size=0.2, random_state=seed)
    train_dataset, test_dataset = split_dataset_by_celltype(dataset, test_size=0.2, random_state=seed)
    print(f"Train samples: {len(train_dataset)}, Test samples: {len(test_dataset)}")
    
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True, num_workers=8)
    test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False, num_workers=8)

    save_path_model_barcode = f"/home/dawn/lab_projects/splicing/for_kai/no_mutagenesis_data/model_barcode_PSI_no_mut_across_cell_types_5epoch_seed{seed}.pth"
    pretrain_barcode_model(device, train_loader, test_loader, epochs, lr, save_path_model_barcode)
    save_path_model_gene_barcode = f"/home/dawn/lab_projects/splicing/for_kai/no_mutagenesis_data/model_barcode_gene_PSI_no_mut_across_cell_types_5epoch_seed{seed}.pth"
    finetune_with_gene_expression(device, train_loader, test_loader, epochs, lr, save_path_model_barcode, save_path_model_gene_barcode)




# ## Evaluation
# model_gene_barcode.load_state_dict(torch.load("DC_TMP_model_barcode_gene_PSI_5epoch_seed0.pth"))
# model_gene_barcode.eval()
# predict_psi = []
# real_psi = []
# psi_residuals = []
# celltypes = []
# barcode_names = []
# with torch.no_grad():
#     psi_loss = 0
#     recon_loss = 0
#     for i, (celltype, barcode_name, barcode, gx, psi) in enumerate(test_loader):
#         barcode = barcode.to(device).float()
#         barcode = barcode.permute(0, 2, 1)
#         gx = gx.to(device).float()
#         psi = psi.to(device).float()
    
#         condition = model_barcode.Barcode_model(barcode)
#         psi_pretrain = model_barcode(barcode)
    
#         psi_residual = model_gene_barcode(barcode, gx, condition)
#         psi_loss += MSELoss(psi_residual.squeeze(), psi-psi_pretrain.squeeze()).item()
        
#         predict_psi.extend(psi_residual.squeeze().cpu().numpy()+psi_pretrain.squeeze().cpu().numpy())
#         psi_residuals.extend(psi_residual.squeeze().cpu().numpy())
#         real_psi.extend(psi.cpu().numpy())
#         celltypes.extend(celltype)
#         barcode_names.extend(barcode_name)


# # In[15]:


# # plot correlation
# import matplotlib.pyplot as plt
# import seaborn as sns

# from scipy.stats import pearsonr, spearmanr
# correlation, _ = pearsonr(real_psi, predict_psi)
# print(f"Pearson Correlation: {correlation:.4f}")
# correlation, _ = spearmanr(real_psi, predict_psi)
# print(f"Spearman Correlation: {correlation:.4f}")
# np.corrcoef(real_psi, predict_psi)[0,1]

# plt.figure(figsize=(12, 12))
# plt.scatter(real_psi, predict_psi, alpha=0.5, s=1)
# plt.xlabel("Real ratio")
# plt.ylabel("Predicted ratio")
# # hide the right and top spines
# plt.gca().spines['right'].set_visible(False)
# plt.gca().spines['top'].set_visible(False)
# plt.show()

