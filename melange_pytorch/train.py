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
from melange_pytorch import Melange

import os
import numpy as np


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


def train(model, train_loader, test_loader, epochs, batch_size, criterion, optimizer):
    
    test_loss_list = []
    for epoch in range(epochs):
        model.eval()
        test_loss = []
        Y_psi_list = []
        output_psi_list = []
        for batch in test_loader:
            X = batch[0].float().cuda()
            Y = batch[1].cuda()
            output = model(X)
            output_psi = output.squeeze()
            Y_psi = Y[:, 1] / (Y[:, 0] + Y[:, 1])
            loss = criterion(output_psi, Y_psi)
            test_loss.append(loss.item())
        print("Epoch: ", epoch, "Test loss: ", np.mean(test_loss))
        test_loss_list.append(np.mean(test_loss))
        torch.cuda.empty_cache()
        
        model.train()
        for i, batch in enumerate(train_loader):
            
            X = batch[0].float().cuda()
            Y = batch[1].cuda()

            optimizer.zero_grad()
            output = model(X)
            # gene_latent = GeneExp(output)
            # output = torch.cat((gene_latent, output), dim=1)
            # latent = Merge(output)
            # The order is [skipped_count, included_count]. PSI is % inclusion.
            output_psi = output.squeeze()
            Y_psi = Y[:, 1] / (Y[:, 0] + Y[:, 1])

            loss = criterion(output_psi, Y_psi)
            loss.backward()
            optimizer.step()
            if i % 100 == 0:
                print("epoch: ", epoch, "batch: ", i, "loss: ", loss.item())

    return test_loss_list

def test(model, test_loader, batch_size):
    model.eval()
    output_psi_list = []
    for batch in test_loader:
        X = batch[0].float().cuda()
        Y = batch[1].cuda()
        output = model(X)
        output_psi = output.squeeze()
        Y_psi = Y[:, 1] / (Y[:, 0] + Y[:, 1])   
        output_psi_list.append(output_psi.cpu().detach())
    output_psi = torch.cat(output_psi_list, dim=0)
    return output_psi

def process_h5py(data):
    shards_idx = np.arange(len(data.keys()) // 2)
    X_list = []
    for idx in shards_idx:
        X = data[f"X{idx}"]
        X_list.append(X)
    X = torch.from_numpy(np.concatenate(X_list, axis=0))
    Y_list = []
    for idx in shards_idx:
        Y = data[f"Y{idx}"]
        Y = Y[0, ...].squeeze()
        Y_list.append(Y)
    Y = torch.from_numpy(np.concatenate(Y_list, axis=0))
    X = X.transpose(1, 2)
    return X, Y
    
if __name__ == "__main__":
    args = parse_args()
    set_seed_for_all(args.seed)
    # Disable wandb if not using it.
    if not args.use_wandb:
        os.environ["WANDB_MODE"] = "disabled"

    # Initialize wandb
    wandb.init(
        project="testing",
        config={
            "learning_rate": args.lr,
            "architecture": args.model,
            "dataset": "spliceai",
            "epochs": args.epochs,
        },
        reinit=True,
    )

    train_data = h5py.File(args.train_h5, "r")
    test_data = h5py.File(args.test_h5, "r")
    
    # Divide by 2 because the dataset has both X and Y. (e.g. X0, Y0, X1, Y1, ... Xn, Yn)
    # So it really just have half the number of sample batches.
    X_train, Y_train = process_h5py(train_data)
    X_test, Y_test = process_h5py(test_data)
    
    train_data = TensorDataset(X_train, Y_train)
    test_data = TensorDataset(X_test, Y_test)    
    
    train_loader = DataLoader(train_data, batch_size=args.batch_size, shuffle=True)
    test_loader = DataLoader(test_data, batch_size=args.batch_size, shuffle=False)

    print("Train data shape: ", X_train.shape)
    print("Test data shape: ", X_test.shape)
    
    # # # Initialize model.
    model = Melange()
    # GeneExp = GeneExp()
    # Merge = Merge()
    model.cuda()

    # Criterion and optimizer.
    # criterion = nn.BCEWithLogitsLoss()
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=1e-4, weight_decay=args.weight_decay)

    test_loss = train(model, train_loader, test_loader, args.epochs, args.batch_size, criterion, optimizer)
    torch.cuda.empty_cache()
    # model.load_state_dict(torch.load(args.model_path))
    output_psi = test(model, test_loader, args.batch_size)
    Y_psi = Y_test[:, 1] / (Y_test[:, 0] + Y_test[:, 1])
    print("Output: ", output_psi[0:10])
    print("Y: ", Y_psi[0:10])
    print("Test loss: ", test_loss[0:10])
    np.savetxt("output_psi.txt", output_psi.cpu().detach().numpy())
    np.savetxt("Y_psi.txt", Y_psi.cpu().detach().numpy())
    np.savetxt("test_loss.txt", test_loss)
    torch.save(model.state_dict(), args.model_path)