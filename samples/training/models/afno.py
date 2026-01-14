import torch
import torch.nn as nn
import torch.nn.functional as F
from functools import partial
from typing import List

torch.manual_seed(943442)
Tensor = torch.Tensor
import os

class Normalizer(torch.nn.Module):
    def __init__(self, size, name, max_accumulation=10**8, std_epsilon=1e-8):
        super(Normalizer, self).__init__()
        self.name = name
        self.max_accumulation = max_accumulation
        self.std_epsilon = std_epsilon

        self.register_buffer('acc_count', torch.zeros(1, dtype=torch.float32))
        self.register_buffer('num_acc', torch.zeros(1, dtype=torch.float32))
        self.register_buffer('acc_sum', torch.zeros(size, dtype=torch.float32))
        self.register_buffer('acc_sum_squared', torch.zeros(size, dtype=torch.float32))


    def forward(self, batched_data: torch.Tensor, accumulate: bool) -> torch.Tensor:
        """Normalizes input data and accumulates statistics."""
        if accumulate and self.num_acc < self.max_accumulation:
            self._accumulate(batched_data)
        batched_data = (batched_data - self._mean()) / self._std_with_eps()
        return batched_data

    def inverse(self, normalized_batch_data):
        """Inverse transformation of the normalizer."""
        return normalized_batch_data * self._std_with_eps() + self._mean()

    def _accumulate(self, batched_data):
        data_sum = torch.sum(batched_data, dim=[2,3,4], keepdim=True)
        data_sum_squared = torch.sum(batched_data**2, dim=[2,3,4], keepdim=True)

        self.acc_sum += data_sum
        self.acc_sum_squared += data_sum_squared
        self.acc_count += torch.tensor(batched_data[0,0,:,:,:].numel()).to(self.num_acc.device)
        self.num_acc += torch.tensor(1.0).to(self.num_acc.device)

    def _mean(self):
        safe_count = torch.maximum(self.acc_count, torch.tensor(1.0).to(self.acc_count.device))
        return self.acc_sum / safe_count

    def _std_with_eps(self):
        safe_count = torch.maximum(self.acc_count, torch.tensor(1.0).to(self.acc_count.device))
        std = torch.sqrt(self.acc_sum_squared / safe_count - self._mean()**2)
        return torch.maximum(std, torch.tensor(self.std_epsilon).to(self.acc_count.device))

class AFNOMlp(nn.Module):
    """Fully-connected Multi-layer perception used inside AFNO

    Parameters
    ----------
    in_features : int
        Input feature size
    latent_features : int
        Latent feature size
    out_features : int
        Output feature size
    activation_fn :  nn.Module, optional
        Activation function, by default nn.GELU
    drop : float, optional
        Drop out rate, by default 0.0
    """

    def __init__(
        self,
        in_features: int,
        latent_features: int,
        out_features: int,
        activation_fn: nn.Module = nn.GELU(),
        drop: float = 0.0,
    ):
        super().__init__()
        self.fc1 = nn.Linear(in_features, latent_features)
        self.act = activation_fn
        self.fc2 = nn.Linear(latent_features, out_features)
        self.drop = nn.Dropout(drop)

    def forward(self, x: Tensor) -> Tensor:
        x = self.fc1(x)
        x = self.act(x)
        x = self.drop(x)
        x = self.fc2(x)
        x = self.drop(x)
        return x

class AFNO3DLayer(nn.Module):
    def __init__(
            self, 
            hidden_size: int, 
            num_blocks: int = 8, 
            sparsity_threshold: float = 0.01,
            hard_thresholding_fraction: float = 1, 
            hidden_size_factor: int =1
        ):

            super().__init__()
            assert hidden_size % num_blocks == 0


            self.hidden_size = hidden_size
            self.sparsity_threshold = sparsity_threshold
            self.num_blocks = num_blocks
            self.block_size = hidden_size // num_blocks
            self.hard_thresholding_fraction = hard_thresholding_fraction
            self.hidden_size_factor = hidden_size_factor
            self.scale = 0.02

            
            self.w1 = nn.Parameter(
                        self.scale 
                        * torch.randn(
                            2, 
                            self.num_blocks, 
                            self.block_size, 
                            self.block_size * hidden_size_factor
                        )
                    )
        

            self.b1 = nn.Parameter(
                        self.scale 
                        * torch.randn(
                            2, 
                            self.num_blocks, 
                            self.block_size * hidden_size_factor
                        )
                    )
            

            self.w2 = nn.Parameter(
                        self.scale 
                        * torch.randn(
                            2, 
                            self.num_blocks, 
                            self.block_size * hidden_size_factor, 
                            self.block_size
                        )
                    )
            
            self.b2 = nn.Parameter(self.scale * torch.randn(2, self.num_blocks, self.block_size)
                                   )
            
    
    def forward(self, x):  # x: [B, D, H, W, C]

        bias = x
        dtype = x.dtype
        x = x.float()
        B, D, H, W, C = x.shape
        dev = x.device

        # compute the fft now

        x = torch.fft.rfftn(x, dim=(1, 2, 3), norm='ortho') 
        x_real = x.real.view(B, D, H, W//2 + 1, self.num_blocks, self.block_size).to(x.device)
        x_imag = x.imag.view(B, D, H, W//2 + 1, self.num_blocks, self.block_size).to(x.device)
        self.w1 = self.w1.to(x.device)
        self.b1 = self.b1.to(x.device)
        self.w2 = self.w2.to(x.device)
        self.b2 = self.b2.to(x.device)
        


        o1_real = torch.zeros(
            [
                B,
                D,
                H,
                W // 2 + 1,
                self.num_blocks,
                self.block_size * self.hidden_size_factor,
            ],
            device=x.device,
        )

        o1_imag = torch.zeros(
            [
                B,
                D,
                H,
                W // 2 + 1,
                self.num_blocks,
                self.block_size * self.hidden_size_factor,
            ],
            device=x.device,
        )

        kept_d = int(D * self.hard_thresholding_fraction)
        kept_h = int(H * self.hard_thresholding_fraction)
        kept_w = int((W // 2 + 1) * self.hard_thresholding_fraction)

        
        o1_real[:, :kept_d, :kept_h, :kept_w] = F.relu(
            torch.einsum("nzyxbi,bio->nzyxbo", x_real[:, :kept_d, :kept_h, :kept_w], self.w1[0]) -
            torch.einsum("nzyxbi,bio->nzyxbo", x_imag[:, :kept_d, :kept_h, :kept_d], self.w1[1]) +
            self.b1[0]
        )

        o1_imag[:, :kept_d, :kept_h, :kept_w] = F.relu(
            torch.einsum("nzyxbi,bio->nzyxbo", x_imag[:, :kept_d, :kept_h, :kept_w], self.w1[0])
            + torch.einsum("nzyxbi,bio->nzyxbo", x_real[:, :kept_d, :kept_h, :kept_w], self.w1[1])
            + self.b1[1]
        )


        o2_real = (
            torch.einsum("nzyxbi,bio->nzyxbo", o1_real, self.w2[0])
            - torch.einsum("nzyxbi,bio->nzyxbo", o1_imag, self.w2[1])
            + self.b2[0]
        )
        o2_imag = (
            torch.einsum("nzyxbi,bio->nzyxbo", o1_imag, self.w2[0])
            + torch.einsum("nzyxbi,bio->nzyxbo", o1_real, self.w2[1])
            + self.b2[1]
        )

        o2 = torch.stack([o2_real[:, :kept_d, :kept_h, :kept_w], o2_imag[:, :kept_d, :kept_h, :kept_w]], dim=-1).to(x.device)
        x = F.softshrink(o2, lambd=self.sparsity_threshold)
        x = torch.view_as_complex(x)
        
        x = x.reshape(B, D, H, W // 2 + 1, C)
        x = torch.fft.irfftn(x, s=(D, H, W), dim=(1, 2, 3), norm="ortho")

        x = x.type(dtype)

        return x.to(dev) + bias.to(dev)

class Block(nn.Module):
    """AFNO block, spectral convolution and MLP

    Parameters
    ----------
    embed_dim : int
        Embedded feature dimensionality
    num_blocks : int, optional
        Number of blocks used in the block diagonal weight matrix, by default 8
    mlp_ratio : float, optional
        Ratio of MLP latent variable size to input feature size, by default 4.0
    drop : float, optional
        Drop out rate in MLP, by default 0.0
    activation_fn: nn.Module, optional
        Activation function used in MLP, by default nn.GELU
    norm_layer : nn.Module, optional
        Normalization function, by default nn.LayerNorm
    double_skip : bool, optional
        Residual, by default True
    sparsity_threshold : float, optional
        Sparsity threshold (softshrink) of spectral features, by default 0.01
    hard_thresholding_fraction : float, optional
        Threshold for limiting number of modes used [0,1], by default 1
    """

    def __init__(
        self,
        embed_dim: int,
        num_blocks: int = 8,
        mlp_ratio: float = 4.0,
        drop: float = 0.0,
        activation_fn: nn.Module = nn.GELU(),
        norm_layer: nn.Module = nn.LayerNorm,
        double_skip: bool = True,
        sparsity_threshold: float = 0.01,
        hard_thresholding_fraction: float = 1.0,
    ):
        super().__init__()
        self.norm1 = norm_layer(embed_dim)
        self.filter = AFNO3DLayer(
            embed_dim, num_blocks, sparsity_threshold, hard_thresholding_fraction
        )
        # self.drop_path = nn.Identity()
        self.norm2 = norm_layer(embed_dim)
        mlp_latent_dim = int(embed_dim * mlp_ratio)
        self.mlp = AFNOMlp(
            in_features=embed_dim,
            latent_features=mlp_latent_dim,
            out_features=embed_dim,
            activation_fn=activation_fn,
            drop=drop,
        )
        self.double_skip = double_skip

    def forward(self, x: Tensor) -> Tensor:
        dev = x.device
        residual = x
        x = self.norm1(x)
        x = self.filter(x)

        if self.double_skip:
            x = x + residual
            residual = x

        x = self.norm2(x)
        x = self.mlp(x)
        x = x + residual
        return x



class PatchEmbed(nn.Module):
    """Patch embedding layer

    Converts 3D patch into a 1D vector for input to AFNO

    Parameters
    ----------
    inp_shape : List[int]
        Input image dimensions [Depth, height, width]
    in_channels : int
        Number of input channels
    patch_size : List[int], optional
        Size of image patches, by default [1, 1, 1]
    embed_dim : int, optional
        Embedded channel size, by default 256
    """

    def __init__(
        self,
        inp_shape: List[int],
        in_channels: int,
        patch_size: List[int] = [1, 1, 1],
        embed_dim: int = 256,
    ):
        super().__init__()
        if len(inp_shape) != 3:
            raise ValueError("inp_shape should be a list of length 3")
        if len(patch_size) != 3:
            raise ValueError("patch_size should be a list of length 3")

        num_patches = (inp_shape[2] // patch_size[2]) * (inp_shape[1] // patch_size[1]) * (inp_shape[0] // patch_size[0])

        self.inp_shape = inp_shape
        self.patch_size = patch_size
        self.num_patches = num_patches

        self.proj = nn.Conv3d(
            in_channels, embed_dim, kernel_size=patch_size, stride=patch_size
        )

    def forward(self, x: Tensor) -> Tensor:
        B, C, D, H, W = x.shape
        """
        if [D, H, W] != self.inp_shape:
            raise ValueError(
                f"Input volume size ({D}*{H}*{W}) doesn't match model ({self.inp_shape})."
            )
        """
        x = self.proj(x).flatten(2).transpose(1, 2)
        return x


class AFNO(nn.Module):
    """Adaptive Fourier neural operator (AFNO) model.

    Note
    ----
    AFNO is a model that is designed for 3D images only.

    Parameters
    ----------
    inp_shape : List[int]
        Input image dimensions [depth, height, width]
    in_channels : int
        Number of input channels
    out_channels: int
        Number of output channels
    patch_size : List[int], optional
        Size of image patches, by default [1, 1, 1]
    embed_dim : int, optional
        Embedded channel size, by default 256
    depth : int, optional
        Number of AFNO layers, by default 4
    mlp_ratio : float, optional
        Ratio of layer MLP latent variable size to input feature size, by default 4.0
    drop_rate : float, optional
        Drop out rate in layer MLPs, by default 0.0
    num_blocks : int, optional
        Number of blocks in the block-diag frequency weight matrices, by default 16
    sparsity_threshold : float, optional
        Sparsity threshold (softshrink) of spectral features, by default 0.01
    hard_thresholding_fraction : float, optional
        Threshold for limiting number of modes used [0,1], by default 1

    Example
    -------
    >>> model = physicsnemo.models.afno.AFNO(
    ...     inp_shape=[38, 38, 38],
    ...     in_channels=3,
    ...     out_channels=6,
    ...     patch_size=(2, 2, 2),
    ...     embed_dim=256,
    ...     depth=4,
    ...     num_blocks=16,
    ... )
    >>> input = torch.randn(N, 3, 38, 38, 38) #(N, C, D, H, W)
    >>> output = model(input)
    >>> output.size()
    torch.Size([N, 6, 38, 38, 38])


    """

    def __init__(
        self,
        inp_shape: List[int] = [38, 38, 38],
        in_channels: int = 3,
        out_channels: int = 6,
        patch_size: List[int] = [1, 1, 1],
        embed_dim: int = 256,
        depth: int = 4,
        mlp_ratio: float = 4.0,
        drop_rate: float = 0.0,
        num_blocks: int = 16,
        sparsity_threshold: float = 0.01,
        hard_thresholding_fraction: float = 1.0,
    ) -> None:
        super().__init__()
        if len(inp_shape) != 3:
            raise ValueError("inp_shape should be a list of length 3")
        if len(patch_size) != 3:
            raise ValueError("patch_size should be a list of length 3")

        if not (
            inp_shape[0] % patch_size[0] == 0 and inp_shape[1] % patch_size[1] == 0
        ):
            raise ValueError(
                f"input shape {inp_shape} should be divisible by patch_size {patch_size}"
            )
        
        self.in_chans = in_channels
        self.out_chans = out_channels
        self.inp_shape = inp_shape
        self.patch_size = patch_size
        self.embed_dim = embed_dim
        self.num_blocks = num_blocks
        norm_layer = partial(nn.LayerNorm, eps=1e-6)
        self.D, self.H, self.W = self.inp_shape
        self.pD, self.pH, self.pW = self.patch_size

        self.normalizer = Normalizer(size=[1, in_channels, 1, 1, 1], name='input')

        
        self.patch_embed = PatchEmbed(
            inp_shape=inp_shape,
            in_channels=self.in_chans,
            patch_size=self.patch_size,
            embed_dim=embed_dim,
        )
        num_patches = self.patch_embed.num_patches

        self.pos_embed = nn.Parameter(torch.zeros(1, num_patches, embed_dim))
        self.pos_drop = nn.Dropout(p=drop_rate)

        self.d = inp_shape[0] // self.patch_size[0]
        self.h = inp_shape[1] // self.patch_size[1]
        self.w = inp_shape[2] // self.patch_size[2]

        self.blocks = nn.ModuleList(
            [
                Block(
                    embed_dim=embed_dim,
                    num_blocks=self.num_blocks,
                    mlp_ratio=mlp_ratio,
                    drop=drop_rate,
                    norm_layer=norm_layer,
                    sparsity_threshold=sparsity_threshold,
                    hard_thresholding_fraction=hard_thresholding_fraction,
                )
                for i in range(depth)
            ]
        )

        self.head = nn.Linear(
            embed_dim,
            self.out_chans * self.patch_size[0] * self.patch_size[1] * self.patch_size[2],
            bias=False,
        )


        torch.nn.init.trunc_normal_(self.pos_embed, std=0.02)
        self.apply(self._init_weights)
        
    def _init_weights(self, m):
        """Init model weights"""
        if isinstance(m, nn.Linear):
            torch.nn.init.trunc_normal_(m.weight, std=0.02)
            if m.bias is not None:
                nn.init.constant_(m.bias, 0)
        elif isinstance(m, nn.LayerNorm):
            nn.init.constant_(m.bias, 0)
            nn.init.constant_(m.weight, 1.0)

    def forward_features(self, x: Tensor) -> Tensor:
        """Forward pass of core AFNO"""
        dev = x.device
        B = x.shape[0]
        x = self.patch_embed(x)
        x = x + self.pos_embed
        x = self.pos_drop(x)

        x = x.reshape(B, self.d, self.h, self.w, self.embed_dim)
        for blk in self.blocks:
            x = blk(x)
        return x
    
    def forward(self, x: Tensor) -> Tensor:
        x = self.normalizer(x, accumulate=True).float() 
        x = self.forward_features(x)
        x = self.head(x)

        # Correct tensor shape back into [B, D, C, H, W]
        # [b d h w (p1 p2 p3 c_out)]
        out = x.view(list(x.shape[:-1]) + [self.patch_size[0], self.patch_size[1], self.patch_size[2], -1])
        #out = x.view(x.shape[0], self.d, self.h, self.w, self.patch_size[0], self.patch_size[1], self.patch_size[2], self.out_chans)
        # [b d h w p1 p2 p3 c_out]
        #out = torch.permute(out, (0, 7, 1, 6, 2, 4, 3, 5))
        out = out.permute(0, 7, 1, 4, 2, 5, 3, 6)
        # [b c_out, d, p3, h, p1, w, p2]
        out = out.reshape(x.size(0), self.out_chans, self.D, self.H, self.W)
        # [b c_out, (h*p1), (w*p2)]
        return out.double()
