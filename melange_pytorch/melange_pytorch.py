import torch
import torch.nn as nn
from einops import rearrange
import torch.nn.functional as F

class Residual(nn.Module):
    def __init__(self,fn):
        super().__init__()
        self.fn = fn
    
    def forward(self, x):
        return x + self.fn(x)
    
def ResidualBlock(in_channels, out_channels, kernel_size, dilation):
    return Residual(nn.Sequential(
        nn.BatchNorm1d(in_channels),
        nn.ReLU(),
        nn.Conv1d(in_channels, out_channels, kernel_size, padding="same", dilation=dilation),
        nn.BatchNorm1d(out_channels),
        nn.ReLU(),
        nn.Conv1d(out_channels, out_channels, kernel_size, padding="same", dilation=dilation)
    )) 


# class Melange(nn.Module):
#     def __init__(self):
#         super().__init__()
#         self.conv1 = nn.Conv1d(4, 32, 1, padding = "same")
#         self.res_conv1 = nn.Conv1d(32, 32, 1, padding = "same")

#         self.block1 = nn.Sequential(
#             ResidualBlock(32, 32, 11, 1),
#             ResidualBlock(32, 32, 11, 1),
#             ResidualBlock(32, 32, 11, 1),
#             ResidualBlock(32, 32, 11, 1),
#             nn.Conv1d(32, 32, 1, dilation = 1, padding = "same")
#         )
#         self.adaptive_pool = nn.AdaptiveAvgPool1d(1)
#         self.conv_last = nn.Conv1d(32, 2, 1, padding = "same")
    
#     def forward(self, x):
#         # print("X shape: ", x.shape)
#         x = self.conv1(x)
#         # print("After conv1 shape: ", x.shape)
#         detour = self.res_conv1(x)
#         # print("Detour shape: ", detour.shape)
#         x = detour + self.block1(x)
#         x = self.adaptive_pool(x)
#         x = self.conv_last(x)
#         # print("Final X shape: ", x.shape)
#         return rearrange(x, 'b c l -> b l c')

class Melange(nn.Module):
    def __init__(self):
        super().__init__()
        # Convolutional layers
        self.conv1 = nn.Conv1d(in_channels=4, out_channels=32, kernel_size=11, padding='same')
        self.pool1 = nn.MaxPool1d(kernel_size=2, stride=2)

        self.conv2 = nn.Conv1d(in_channels=32, out_channels=32, kernel_size=11, padding='same')
        self.pool2 = nn.MaxPool1d(kernel_size=2, stride=2)

        self.conv3 = nn.Conv1d(in_channels=32, out_channels=2, kernel_size=11, padding='same')
        self.pool3 = nn.MaxPool1d(kernel_size=2, stride=2)

        # Linear layers
        # Assuming the sequence length is properly reduced to 1 after convolutions and pooling
        self.fc1 = nn.Linear(in_features=224, out_features=32)
        # self.fc2 = nn.Linear(in_features=32, out_features=2)  # To get an output of 2 features

    def forward(self, x):
        # Apply convolutional layers with ReLU activations and max pooling
        x = F.relu(self.conv1(x))
        x = self.pool1(x)

        x = F.relu(self.conv2(x))
        x = self.pool2(x)

        x = F.relu(self.conv3(x))
        x = self.pool3(x)

        # Flatten the output for the linear layer
        x = x.view(x.size(0), -1)  # Flatten
        
        # Apply linear layers with ReLU activation
        x = F.relu(self.fc1(x))
        # x = self.fc2(x)

        # Reshape to get the desired output shape (batch, 1, 2)
        x = x.view(-1, 1, 1)

        return x
    
if __name__ == "__main__":
    model = Melange()
    x = torch.randn(32, 4, 900) # Shape is [BATCH_SIZE, 4 (ACGT), SEQUENCE_LENGTH]
    # print(x)
    y = model(x)
    print(y) # This is [1, BATCH_SIZE, SEQUENCE_LENGTH, 2]