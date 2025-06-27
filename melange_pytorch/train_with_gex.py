"""
python train.py --output train_output.torchfile \
    --train_h5 ../data/dataset_train_all.h5 \
        --test_h5 ../data/dataset_test_0.h5 --epochs 1

"""

import argparse
import tqdm
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


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=str, default="SpliceAI")
    parser.add_argument("--output", required=True)
    parser.add_argument("--train_h5", required=True)
    parser.add_argument("--test_h5", required=True)
    parser.add_argument("--use-wandb", action="store_true", default=False)

    parser.add_argument("--epochs", type=int, default=10)
    parser.add_argument("--batch_size", type=int, default=32)

    parser.add_argument("--lr", type=float, default=0.001)
    parser.add_argument("--weight_decay", type=float, default=0.0001)
    parser.add_argument(
        "--model_path", type=str, default="melange_pytorch/model/model.pt"
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


def train(
    sequence_cnn_model,
    gex_model,
    merge_model,
    h5f,
    gex_df,
    train_idxes,
    batch_size,
    criterion,
    sequence_cnn_optimizer,
    gex_optimizer,
    merge_optimizer,
):
    sequence_cnn_model.train()
    gex_model.train()
    merge_model.train()
    running_output, running_label = [], []

    batch_idx = 0
    print("Number of train idxes:", len(train_idxes))
    for i, train_idx in enumerate(
        np.random.choice(train_idxes, size=len(train_idxes), replace=False)
    ):
        # if i > 1:
        #     break
        batch_idx += 1
        # The initial shape is (batch_size, 900, 4)
        X = h5f[f"X{train_idx}"]
        # Need to rearrange to (batch_size, 4, 900)
        X = rearrange(X[:], "b c l -> b l c")
        # This has shape (1, batch_size, 900, 2)
        Y = h5f[f"Y{train_idx}"]
        # Just keep the first dimension.
        # Final shape is (batch_size, 900, 2)
        Y = Y[0, ...]

        Z = h5f[f"Z{train_idx}"][:]
        Z_array = []
        for j in range(len(Z)):
            celltype = Z[j].decode("utf-8")
            # print(celltype)
            gex = gex_df.loc[celltype].values.copy()
            # Reshape this to a 2D array.
            gex = gex.reshape(1, -1)
            # print(gex.shape)
            Z_array.extend([gex])
        Z_array = np.array(Z_array)

        dataset = TensorDataset(
            torch.from_numpy(X).float(),
            torch.from_numpy(Y).float(),
            torch.from_numpy(Z_array).float(),
        )
        loader = DataLoader(
            dataset,
            batch_size=batch_size,
            shuffle=True,
            num_workers=8,
            pin_memory=True,
            drop_last=True,
        )

        bar = tqdm.tqdm(loader, desc=f"Shard {i}/{len(train_idxes)}", leave=False)

        for batch in bar:
            X = batch[0].cuda()
            Y = batch[1].cuda()
            Z = batch[2].cuda()

            sequence_cnn_optimizer.zero_grad()
            gex_optimizer.zero_grad()
            merge_optimizer.zero_grad()

            sequence_output = sequence_cnn_model(X)
            gene_latent = gex_model(Z)
            output = torch.cat((gene_latent, sequence_output), dim=2)
            output = merge_model(output)
            # The order is [skipped_count, included_count]. PSI is % inclusion.
            # ouput_psi = output[:, :, 1] / (output[:, :, 0] + output[:, :, 1])
            output_psi = output
            # print(ouput_psi.cpu().detach().numpy())
            Y_psi = Y[:, :, 1] / (Y[:, :, 0] + Y[:, :, 1])

            loss = criterion(output_psi, Y_psi)
            loss.backward()
            sequence_cnn_optimizer.step()
            gex_optimizer.step()
            merge_optimizer.step()

            running_output.append(output.detach().cpu())
            running_label.append(Y.detach().cpu())

            if batch_idx % 20 == 0:
                running_output = torch.cat(running_output, dim=0)
                running_label = torch.cat(running_label, dim=0)
                running_output_psi = running_output
                running_label_psi = running_label[:, :, 1] / (
                    running_label[:, :, 0] + running_label[:, :, 1]
                )

                running_loss = criterion(running_output_psi, running_label_psi)
                wandb.log(
                    {
                        "train_loss": running_loss.item(),
                    }
                )

                running_output, running_label = [], []


def test(
    sequence_cnn_model,
    gex_model,
    merge_model,
    h5f,
    gex_df,
    test_shard_ids,
    batch_size,
    criterion,
):
    sequence_cnn_model.eval()
    gex_model.eval()
    merge_model.eval()
    out, label = [], []
    for shard_idx in test_shard_ids:
        X = h5f[f"X{shard_idx}"]
        # Need to rearrange to (batch_size, 4, 900)
        X = rearrange(X[:], "b c l -> b l c")
        # This has shape (1, batch_size, 900, 2)
        Y = h5f[f"Y{shard_idx}"]
        # Just keep the first dimension.
        # Final shape is (batch_size, 900, 2)
        Y = Y[0, ...]
        
        Z = h5f[f"Z{shard_idx}"][:]
        Z_array = []
        for i in range(len(Z)):
            celltype = Z[i].decode("utf-8")
            # print(celltype)
            gex = gex_df.loc[celltype].values.copy()
            # Reshape this to a 2D array.
            gex = gex.reshape(1, -1)
            # print(gex.shape)
            Z_array.extend([gex])
        Z_array = np.array(Z_array)

        dataset = TensorDataset(
            torch.from_numpy(X).float(), torch.from_numpy(Y).float(), torch.from_numpy(Z_array).float()
        )
        loader = DataLoader(
            dataset,
            batch_size=batch_size,
            shuffle=False,
            num_workers=8,
            pin_memory=True,
        )
        bar = tqdm.tqdm(loader, leave=False, total=len(loader))
        for batch in bar:
            X = batch[0].cuda()
            Y = batch[1].cuda()
            Z = batch[2].cuda()
            sequence_output = sequence_cnn_model(X)
            gene_latent = gex_model(Z)
            output = torch.cat((gene_latent, sequence_output), dim=2)
            _out = merge_model(output).detach().cpu()
            _label = Y[:, :, 1] / (Y[:, :, 0] + Y[:, :, 1])
            _label = _label.detach().cpu()
            out.append(_out)
            label.append(_label)
    out = torch.cat(out, dim=0)
    label = torch.cat(label, dim=0)
    loss = criterion(out, label)
    print("Test loss: ", loss.item())
    wandb.log({"test_loss": loss.item()})
    return loss.item()


def validate(
    sequence_cnn_model,
    gex_model,
    merge_model,
    h5f,
    gex_df,
    val_shard_ids,
    batch_size,
    criterion,
):
    sequence_cnn_model.eval()
    gex_model.eval()
    merge_model.eval()
    out, label = [], []
    for shard_idx in val_shard_ids:
        X = h5f[f"X{shard_idx}"]
        # Need to rearrange to (batch_size, 4, 900)
        X = rearrange(X[:], "b c l -> b l c")
        # This has shape (1, batch_size, 900, 2)
        Y = h5f[f"Y{shard_idx}"]
        # Just keep the first dimension.
        # Final shape is (batch_size, 900, 2)
        Y = Y[0, ...]
        
        Z = h5f[f"Z{shard_idx}"][:]
        Z_array = []
        for i in range(len(Z)):
            celltype = Z[i].decode("utf-8")
            # print(celltype)
            gex = gex_df.loc[celltype].values.copy()
            # Reshape this to a 2D array.
            gex = gex.reshape(1, -1)
            # print(gex.shape)
            Z_array.extend([gex])
        Z_array = np.array(Z_array)

        dataset = TensorDataset(
            torch.from_numpy(X).float(), torch.from_numpy(Y).float(), torch.from_numpy(Z_array).float()
        )
        loader = DataLoader(
            dataset,
            batch_size=batch_size,
            shuffle=False,
            num_workers=8,
            pin_memory=True,
        )
        bar = tqdm.tqdm(loader, leave=False, total=len(loader))
        for batch in bar:
            X = batch[0].cuda()
            Y = batch[1].cuda()
            Z = batch[2].cuda()
            sequence_output = sequence_cnn_model(X)
            gene_latent = gex_model(Z)
            output = torch.cat((gene_latent, sequence_output), dim=2)
            _out = merge_model(output).detach().cpu()
            _label = Y[:, :, 1] / (Y[:, :, 0] + Y[:, :, 1])
            _label = _label.detach().cpu()
            out.append(_out)
            label.append(_label)
    out = torch.cat(out, dim=0)
    label = torch.cat(label, dim=0)
    loss = criterion(out, label)
    print("Validation loss: ", loss.item())
    wandb.log({"val_loss": loss.item()})
    return loss.item()


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

    # Also read in the gene expression file.
    ccle_gex_file = "/home/jupyter/melange/data/ccle_gex_shortlist_celllines.csv"
    ccle_gex_file = pd.read_csv(ccle_gex_file)
    # Make the column "StrippedName" the index and drop the column.
    ccle_gex_file = ccle_gex_file.set_index("StrippedName")

    # Divide by 3 because the dataset has both X and Y. (e.g. X0, Y0, X1, Y1, ... Xn, Yn)
    # So it really just have half the number of sample batches.
    num_shards = len(train_data.keys()) // 3
    shards_idx = np.random.permutation(num_shards)
    train_idx = shards_idx[: int(num_shards * 0.9)]
    val_idx = shards_idx[int(num_shards * 0.9) :]

    test_shards_idx = np.arange(len(test_data.keys()) // 2)

    # Initialize model.
    sequence_cnn_model = SequenceCNN()
    gex_model = GexFullyConnected()
    merge_model = MergSequenceAndGex()

    sequence_cnn_model.cuda()
    gex_model.cuda()
    merge_model.cuda()

    # Criterion and optimizer.
    # criterion = nn.CrossEntropyLoss()
    criterion = nn.MSELoss()
    sequence_cnn_optimizer = optim.Adam(sequence_cnn_model.parameters(), lr=1e-3)
    gex_optimizer = optim.Adam(gex_model.parameters(), lr=1e-3)
    merge_optimizer = optim.Adam(merge_model.parameters(), lr=1e-3)

    sequence_cnn_scheduler = optim.lr_scheduler.MultiStepLR(
        sequence_cnn_optimizer, milestones=[6, 7, 8, 9], gamma=0.5
    )
    gex_scheduler = optim.lr_scheduler.MultiStepLR(
        gex_optimizer, milestones=[6, 7, 8, 9], gamma=0.5
    )
    merge_scheduler = optim.lr_scheduler.MultiStepLR(
        merge_optimizer, milestones=[6, 7, 8, 9], gamma=0.5
    )

    for epoch in range(args.epochs):
        train(
            sequence_cnn_model,
            gex_model,
            merge_model,
            train_data,
            ccle_gex_file,
            train_idx,
            args.batch_size,
            criterion,
            sequence_cnn_optimizer,
            gex_optimizer,
            merge_optimizer,
        )
        validate(
            sequence_cnn_model,
            gex_model,
            merge_model,
            train_data,
            ccle_gex_file,
            val_idx,
            args.batch_size,
            criterion,
        )
        test(
            sequence_cnn_model,
            gex_model,
            merge_model,
            test_data,
            ccle_gex_file,
            test_shards_idx,
            args.batch_size,
            criterion,
        )
        
        sequence_cnn_scheduler.step()
        gex_scheduler.step()
        merge_scheduler.step()

    # Save the output.
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    torch.save(sequence_cnn_model.state_dict(), os.path.join(args.output, "sequence_cnn_model.pt"))
    torch.save(gex_model.state_dict(), os.path.join(args.output, "gex_model.pt"))
    torch.save(merge_model.state_dict(), os.path.join(args.output, "merge_model.pt"))
