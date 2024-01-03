import torch
import torch.nn as nn
from einops import rearrange

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


class Melange(nn.Module):
    def __init__(self):
        super().__init__()
        self.conv1 = nn.Conv1d(4, 32, 1, padding = "same")
        self.res_conv1 = nn.Conv1d(32, 32, 1, padding = "same")

        self.block1 = nn.Sequential(
            ResidualBlock(32, 32, 11, 1),
            ResidualBlock(32, 32, 11, 1),
            ResidualBlock(32, 32, 11, 1),
            ResidualBlock(32, 32, 11, 1),
            nn.Conv1d(32, 32, 1, dilation = 1, padding = "same")
        )
        self.adaptive_pool = nn.AdaptiveAvgPool1d(1)
        self.conv_last = nn.Conv1d(32, 2, 1, padding = "same")
    
    def forward(self, x):
        # print("X shape: ", x.shape)
        x = self.conv1(x)
        # print("After conv1 shape: ", x.shape)
        detour = self.res_conv1(x)
        # print("Detour shape: ", detour.shape)
        x = detour + self.block1(x)
        x = self.adaptive_pool(x)
        x = self.conv_last(x)
        # print("Final X shape: ", x.shape)
        return rearrange(x, 'b c l -> b l c')

if __name__ == "__main__":
    model = Melange()
    x = torch.randn(32, 4, 900)
    y = model(x)
    print(y.shape) # This is [BATCH_SIZE, SEQUENCE_LENGTH, 2]