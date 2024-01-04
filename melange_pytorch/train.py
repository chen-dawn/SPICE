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


def train(model, h5f, train_idxes, batch_size, criterion, optimizer):
    model.train()
    running_output, running_label = [], []

    batch_idx = 0
    print("Number of train idxes:", len(train_idxes))
    for i, train_idx in enumerate(
        np.random.choice(train_idxes, size=len(train_idxes), replace=False)
    ):
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

        dataset = TensorDataset(
            torch.from_numpy(X).float(), torch.from_numpy(Y).float()
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

            optimizer.zero_grad()
            output = model(X)
            print(output)
            # The order is [skipped_count, included_count]. PSI is % inclusion.
            # ouput_psi = output[:, :, 1] / (output[:, :, 0] + output[:, :, 1])
            output_psi = output
            Y_psi = Y[:, :, 1] / (Y[:, :, 0] + Y[:, :, 1])
            # print(ouput_psi.cpu().detach().numpy())
            loss = criterion(output_psi, Y_psi)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=5)
            optimizer.step()

            running_output.append(output.detach().cpu())
            running_label.append(Y.detach().cpu())

            if batch_idx % 10 == 0:
                running_output = torch.cat(running_output, dim=0)
                running_label = torch.cat(running_label, dim=0)
                running_output_psi = running_output[:, :, 1] / (
                    running_output[:, :, 0] + running_output[:, :, 1]
                )
                running_label_psi = running_label[:, :, 1] / (
                    running_label[:, :, 0] + running_label[:, :, 1]
                )

                # top_k_acc1 = top_k_accuracy(running_output_prob[:,:,1], running_label[:,:,1])
                # top_k_acc2 = top_k_accuracy(running_output_prob[:,:,2], running_label[:,:,2])
                # auprc_1 = average_precision_score(running_label[:,:,1].view(-1), running_output_prob[:,:,1].view(-1))
                # auprc_2 = average_precision_score(running_label[:,:,2].view(-1), running_output_prob[:,:,2].view(-1))

                running_loss = criterion(running_output_psi, running_label_psi)
                wandb.log(
                    {
                        "train_loss": running_loss.item(),
                        # "top_k_acc1": top_k_acc1.item(),
                        # "top_k_acc2": top_k_acc2.item(),
                        # "auprc_1": auprc_1,
                        # "auprc_2": auprc_2,
                    }
                )
                # wandb.log({
                #     "train_loss": running_loss.item(),
                #     "top_k_acc1": top_k_acc1.item(),
                #     "top_k_acc2": top_k_acc2.item(),
                #     "auprc_1": auprc_1,
                #     "auprc_2": auprc_2,
                # })
                # bar.set_postfix({
                #     "loss": running_loss.item(),
                #     "top_k_acc1": top_k_acc1.item(),
                #     "top_k_acc2": top_k_acc2.item(),
                #     "auprc_1": auprc_1,
                #     "auprc_2": auprc_2,
                # })

                running_output, running_label = [], []


def test(model, h5f, test_shard_ids, batch_size, criterion):
    model.eval()
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


if __name__ == "__main__":
    args = parse_args()
    set_seed_for_all(args.seed)
    if not args.use_wandb:
        os.environ["WANDB_MODE"] = "disabled"

    # Initialize wandb
    wandb.init(
        project="spliceai",
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
    # train_loader = DataLoader(train_data, batch_size=args.batch_size, shuffle=True)
    # test_loader = DataLoader(test_data, batch_size=args.batch_size, shuffle=True)

    # Divide by 2 because the dataset has both X and Y. (e.g. X0, Y0, X1, Y1, ... Xn, Yn)
    # So it really just have half the number of sample batches.
    num_shards = len(train_data.keys()) // 2
    shards_idx = np.random.permutation(num_shards)
    train_idx = shards_idx[: int(num_shards * 0.9)]
    val_idx = shards_idx[int(num_shards * 0.9) :]

    test_shards_idx = np.arange(len(test_data.keys()) // 2)

    # Initialize model.
    model = Melange()
    model.cuda()

    # Criterion and optimizer.
    # criterion = nn.CrossEntropyLoss()
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=1e-3)
    # scheduler = optim.lr_scheduler.MultiStepLR(
    #     optimizer, milestones=[6, 7, 8, 9], gamma=0.5
    # )
    scheduler = optim.lr_scheduler.MultiStepLR(
        optimizer, milestones=[30, 60, 90], gamma=0.1
    )

    for epoch in range(args.epochs):
        train(model, train_data, train_idx, args.batch_size, criterion, optimizer)
        # validate(model, train_data, val_idx, criterion, epoch)
        # test(model, test_data, test_shards_idx, criterion, epoch)

        scheduler.step()

    torch.save(model.state_dict(), args.output)
