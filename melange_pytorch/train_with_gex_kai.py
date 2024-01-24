"""
python train.py --output train_output.torchfile \
    --train_h5 ../data/dataset_train_all.h5 \
        --test_h5 ../data/dataset_test_0.h5 --epochs 1

"""

import argparse
from tqdm import tqdm
import h5py
import random
import wandb  # Library for logging metrics.
from einops import rearrange

import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
from sklearn.metrics import average_precision_score

# Import dataloader
from torch.utils.data import TensorDataset, DataLoader
from melange_pytorch import SequenceCNN, GexFullyConnected, MergSequenceAndGex

import os
import numpy as np

import pandas as pd
import h5py

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=str, default="SpliceAI")
    parser.add_argument("--output", required=True)
    parser.add_argument("--train_h5", required=True)
    parser.add_argument("--test_h5", required=True)
    parser.add_argument("--use-wandb", action="store_true", default=False)

    parser.add_argument("--epochs", type=int, default=100)
    parser.add_argument("--batch_size", type=int, default=256)

    parser.add_argument("--lr", type=float, default=0.001)
    parser.add_argument("--weight_decay", type=float, default=0.0001)
    parser.add_argument(
        "--model_path", type=str, default="model/model.pt"
    )
    parser.add_argument("--log_path", type=str, default="melange_pytorch/log/log.txt")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--device", type=str, default="cuda")

    args = parser.parse_args()
    return args


def shuffle(arr):
    """
    Shuffle the array in-place.
    """
    return np.random.choice(arr, size=len(arr), replace=False)


def top_k_accuracy(pred_probs, labels):
    """
    Predict the top k accuracy. Here the k is not a hyperparameter, but the number of positive labels.
    """
    pred_probs, labels = map(lambda x: x.view(-1), [pred_probs, labels])
    k = (labels == 1.0).sum().item()
    _, top_k_indices = pred_probs.topk(k)
    correct = (labels[top_k_indices] == 1.0).float().mean()
    return correct


def set_seed_for_all(seed):
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    random.seed(seed)

    # Performance drops, so commenting out for now.
    # torch.backends.cudnn.benchmark = False
    # torch.backends.cudnn.deterministic = True


def train(seq_model, gene_model, merge_model, train_loader, test_loader, epochs, batch_size, criterion, optimizer):
    
    test_loss_list = []
    for epoch in range(epochs):
        
        seq_model.train()
        gene_model.train()
        merge_model.train()
        for i, batch in enumerate(train_loader):
            
            X = batch[0].float().cuda()
            Y = batch[1].cuda()
            Z = batch[2].float().cuda()

            optimizer.zero_grad()
            z_seq = seq_model(X)
            z_gene = gene_model(Z)
            z_merge = merge_model(torch.cat((z_seq, z_gene), dim=1))
            output_psi = z_merge.squeeze()
            Y_psi = Y[:, 1] / (Y[:, 0] + Y[:, 1])

            loss = criterion(output_psi, Y_psi)
            loss.backward()
            optimizer.step()
            if i % 100 == 0:
                print("epoch: ", epoch, "batch: ", i, "loss: ", loss.item())
                
        seq_model.eval()
        gene_model.eval()
        merge_model.eval()
        test_loss = []
        Y_psi_list = []
        output_psi_list = []
        for batch in test_loader:
            X = batch[0].float().cuda()
            Y = batch[1].cuda()
            Z = batch[2].float().cuda()
            z_seq = seq_model(X)
            z_gene = gene_model(Z)
            z_merge = merge_model(torch.cat((z_seq, z_gene), dim=1))
            output_psi = z_merge.squeeze()
            Y_psi = Y[:, 1] / (Y[:, 0] + Y[:, 1])
            print(Y_psi)
            loss = criterion(output_psi, Y_psi)
            test_loss.append(loss.item())
        print("Epoch: ", epoch, "Test loss: ", np.mean(test_loss))
        test_loss_list.append(np.mean(test_loss))
        torch.cuda.empty_cache()
    
        torch.save(seq_model.state_dict(), "model/seq_model.pt")
        torch.save(gene_model.state_dict(), "model/gene_model.pt")
        torch.save(merge_model.state_dict(), "model/merge_model.pt")

    return test_loss_list

def test(seq_model, gene_model, merge_model, test_loader, batch_size):
    seq_model.eval()
    gene_model.eval()
    merge_model.eval()
    output_psi_list = []
    for batch in test_loader:
        X = batch[0].float().cuda()
        Y = batch[1].cuda()
        Z = batch[2].float().cuda()
        z_seq = seq_model(X)
        z_gene = gene_model(Z)
        z_merge = merge_model(torch.cat((z_seq, z_gene), dim=1))
        output_psi = z_merge.squeeze()
        Y_psi = Y[:, 1] / (Y[:, 0] + Y[:, 1])   
        output_psi_list.append(output_psi.cpu().detach())
    output_psi = torch.cat(output_psi_list, dim=0)
    return output_psi

def process_h5py(data):
    
    # Also read in the gene expression file.
    ccle_gex_file = "/home/kcao/Dawn/data/ccle_gex_shortlist_celllines.csv"
    ccle_gex_file = pd.read_csv(ccle_gex_file)
    # Make the column "StrippedName" the index and drop the column.
    ccle_gex_file = ccle_gex_file.set_index("StrippedName")
    print(ccle_gex_file)
    
    shards_idx = np.arange(len(data.keys()) // 3)
    # shards_idx = shuffle(shards_idx)
    shards_idx = shards_idx
    
    celltype_to_gene = {}
    celltypes = ['A375', 'GAMG', 'K562', 'KELLY', 'MCF7', 'SKNAS', 'T47D', 'TOV21G', 'U251MG']
    for celltype in celltypes:
        celltype_to_gene[celltype] = ccle_gex_file.loc[celltype].values.copy()
    Z_array = []
    for idx in tqdm(shards_idx):
        Z_list = []
        Z = data[f"Z{idx}"][:]
        for i in range(len(Z)):
            celltype = Z[i].decode("utf-8")
            gex = celltype_to_gene[celltype]
            # Reshape this to a 2D array.
            gex = gex.reshape(1, -1).astype(np.float32)
            # print(gex.shape)
            Z_list.append(gex)
            
        Z_list = np.concatenate(Z_list, axis=0)
        Z_array.append(Z_list)
        
    Z_array = np.concatenate(Z_array, axis=0)
    Z = torch.from_numpy(Z_array)
    
    X_list = []
    for idx in shards_idx:
        X = data[f"X{idx}"]
        X_list.append(X)
    X = torch.from_numpy(np.concatenate(X_list, axis=0))
    X = X.transpose(1, 2)
    
    Y_list = []
    for idx in shards_idx:
        Y = data[f"Y{idx}"]
        Y = Y[0, ...].squeeze()
        Y_list.append(Y)
    Y = torch.from_numpy(np.concatenate(Y_list, axis=0))
        
    return X, Y, Z
    
if __name__ == "__main__":
    args = parse_args()
    set_seed_for_all(args.seed)

    train_data = h5py.File(args.train_h5, "r")
    test_data = h5py.File(args.test_h5, "r")

    # Divide by 3 because the dataset has both X and Y. (e.g. X0, Y0, Z0, X1, Y1, Z1,... Xn, Yn, Zn)
    # So it really just have half the number of sample batches.
    X_train, Y_train, Z_train = process_h5py(train_data)
    X_test, Y_test, Z_test = process_h5py(test_data)
    
    train_data = TensorDataset(X_train, Y_train, Z_train)
    test_data = TensorDataset(X_test, Y_test, Z_test)    
    
    train_loader = DataLoader(train_data, batch_size=args.batch_size, shuffle=True)
    test_loader = DataLoader(test_data, batch_size=args.batch_size, shuffle=False)

    print("Train data shape: ", X_train.shape)
    print("Train Gene data shape: ", Z_train.shape)
    print("Test data shape: ", X_test.shape)
    
    # # # Initialize model.
    seq_model = SequenceCNN().cuda()
    gene_model = GexFullyConnected().cuda()
    merge_model = MergSequenceAndGex().cuda()
    
    seq_model.load_state_dict(torch.load("model/seq_model.pt"))
    gene_model.load_state_dict(torch.load("model/gene_model.pt"))
    merge_model.load_state_dict(torch.load("model/merge_model.pt"))
    
    # Criterion and optimizer.
    criterion = nn.MSELoss()
    #optimize all cnn parameters
    optimizer = optim.Adam(list(seq_model.parameters()) + list(gene_model.parameters()) + list(merge_model.parameters()), lr=args.lr, weight_decay=args.weight_decay)

    test_loss = train(seq_model, gene_model, merge_model, train_loader, test_loader, args.epochs, args.batch_size, criterion, optimizer)
    torch.cuda.empty_cache()
    # model.load_state_dict(torch.load(args.model_path))
    output_psi = test(seq_model, gene_model, merge_model, test_loader, args.batch_size)
    Y_psi = Y_test[:, 1] / (Y_test[:, 0] + Y_test[:, 1])
    print("Output: ", output_psi[0:10])
    print("Y: ", Y_psi[0:10])
    print("Test loss: ", test_loss[0:10])
    np.savetxt("output_psi.txt", output_psi.cpu().detach().numpy())
    np.savetxt("Y_psi.txt", Y_psi.cpu().detach().numpy())
    np.savetxt("test_loss.txt", test_loss)
    
    torch.save(seq_model.state_dict(), "model/seq_model.pt")
    torch.save(gene_model.state_dict(), "model/gene_model.pt")
    torch.save(merge_model.state_dict(), "model/merge_model.pt")
