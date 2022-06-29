from torch import nn
import torch.nn.functional as F

class scEpiLock(nn.Module):
    def __init__(self, n_class):
        super(scEpiLock, self).__init__()
        self.Conv1 = nn.Conv1d(in_channels=4, out_channels=320, kernel_size=8)
        self.Conv2 = nn.Conv1d(in_channels=320, out_channels=480, kernel_size=8)
        self.Conv3 = nn.Conv1d(in_channels=480, out_channels=960, kernel_size=4)
        self.Conv4 = nn.Conv1d(in_channels=960, out_channels=1024, kernel_size=4)
        self.Maxpool = nn.MaxPool1d(kernel_size=4, stride=4)
        self.Drop1 = nn.Dropout(p=0.2)
        self.Drop2 = nn.Dropout(p=0.5)
        self.Linear1 = nn.Linear(57*1024, 925)
        self.Linear2 = nn.Linear(925, n_class)

    def forward(self, input):
        # Convolution Layer 1
        # Input Tensor Shape: [batch_size, 4, 1000]
        # Output Tensor Shape: [batch_size, 320, 993]
        x = self.Conv1(input)
        x = F.relu(x)
        # # Pooling Layer 1
        # # Input Tensor Shape: [batch_size, 320, 993]
        # # Output Tensor Shape: [batch_size, 320, 248]
        # x = self.Maxpool(x)
        # x = self.Drop1(x)

        # Convolution Layer 2
        # Input Tensor Shape: [batch_size, 320, 993]
        # Output Tensor Shape: [batch_size, 480, 986]
        x = self.Conv2(x)
        x = F.relu(x)
        # Pooling Layer 1320
        # Input Tensor Shape: [batch_size, 480, 986]
        # Output Tensor Shape: [batch_size, 480, 246]
        x = self.Maxpool(x)
        x = self.Drop1(x)

        # Convolution Layer 3
        # Input Tensor Shape: [batch_size, 480, 246]
        # Output Tensor Shape: [batch_size, 960, 243]
        x = self.Conv3(x)
        x = F.relu(x)
        # Pooling Layer 2
        # Input Tensor Shape: [batch_size, 480, 243]
        # Output Tensor Shape: [batch_size, 480, 60]
        x = self.Maxpool(x)
        x = self.Drop2(x)

        # Convolution Layer 4
        # Input Tensor Shape: [batch_size, 960, 60]
        # Output Tensor Shape: [batch_size, 1024, 57]
        x = self.Conv4(x)
        x = F.relu(x)
        x = self.Drop2(x)

        # Pooling Layer 3
        # Input Tensor Shape: [batch_size, 1024, 57]
        # Output Tensor Shape: [batch_size, 1024, 14]
        #x = self.Maxpool(x)

        x = x.view(-1, 57*1024)
        x = self.Linear1(x)
        x = F.relu(x)
        x = self.Linear2(x)
        return x