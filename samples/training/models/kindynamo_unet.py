# File implementing an adapted version of U-Net for the kin-dynamo experiments.
# Concept from Ronneberger et al. (2015): "U-Net: Convolutional Networks for Biomedical Image Segmentation".


import torch
import torch.nn as nn
import torch.nn.functional as F


# TODO: Create all missing class and method docstrings
# TODO: Complete all incomplete class and method docstrings (include both class or method descriptions and
#   parameter definitions)


class KinDynamo_UNet(nn.Module):
    """

    # TODO: Update the below input shape and description to match the fact that there are now two separate inputs
    Input Shape: (b, c_in, z, y, x),

    where b is batch size; c_in is number of input channels and is equal to the number of spatial
    dimensions of each mean velocity field vector plus the number of spatial dimensions of each mean
    magnetic field vector; and z, y, and x are the number of grid points of the magnetic vector field
    in the z-, y-, and x-dimensions, respectively. The channels should be in the following order, from
    first to last: magnetic z, magnetic y, magnetic x, velocity z, velocity y, velocity x.

    Output Shape: (b, c_out, z, y, x),

    where b is batch size, c_out is number of output channels and is equal to the number of spatial
    dimensions of each mean residual EMF vector; and z, y, and x are the number of grid points of the
    mean residual EMF vector field in the z-, y-, and x-dimensions, respectively and are equal to the
    number of grid points of the mean velocity vector field (which is also equal to the number of grid
    points of the mean magnetic vector field) in each of these dimensions, respectively. The channels
    should be in the following order, from first to last: mean residual EMF z, mean residual EMF y,
    mean residual EMF x.


    NOTES:

    NOTE: This model is a modified version of U-Net. It is implemented for 3D input instead of
        2D input, and it has a smaller depth (i.e., contains fewer encoder layers and correspondingly
        fewer decoder layers) because the dimensions of our vector field inputs are smaller than the
        dimensions (572 x 572) of U-Net's image inputs. Additionally, to avoid having to crop the
        encoder-side inputs to the skip connection, we perform the upsampling after each deconvolutional
        layer instead of before so that we can easily match the size of the encoder-side inputs on the
        decoder side.

    """
    
    # def __init__(self, in_channels=6, out_channels=3, second_level_channels=16):
    def __init__(self, in_channels, out_channels, second_level_channels):

        super(KinDynamo_UNet, self).__init__()

        channels = second_level_channels

        # ENCODER
        self.encoder_convblock_1 = KinDynamo_UNet._block(in_channels, channels)
        self.downsample_1 = nn.MaxPool3d(kernel_size=2, stride=2)

        self.encoder_convblock_2 = KinDynamo_UNet._block(channels, channels * 2)
        self.downsample_2 = nn.MaxPool3d(kernel_size=2, stride=2)

        # BOTTLENECK
        self.bottleneck = KinDynamo_UNet._block(channels * 2, channels * 4)

        # DECODER
        self.deconv_1 = nn.ConvTranspose3d(channels * 4, channels * 2, kernel_size=2, stride=2)

        # NOTE: Upsampling is implemented in forward() using torch.nn.functional to support model inputs
        #   of varying size.
        self.decoder_convblock_1 = KinDynamo_UNet._block(channels * 4, channels * 2)

        self.deconv_2 = nn.ConvTranspose3d(channels * 2, channels, kernel_size=2, stride=2)
        # NOTE: Upsampling is implemented in forward() using torch.nn.functional to support model inputs
        #   of varying size.
        self.decoder_convblock_2 = KinDynamo_UNet._block(channels * 2, channels)

        self.conv = nn.Conv3d(channels, out_channels, kernel_size=1)


    def forward(self, U_x, B_x):
        """
        Perform a forward pass of the neural network. (I.e., perform inference on the input sample.)

        :param U_x: Mean velocity field values of the input sample
        
        :param B_x: Mean magnetic field values of the input sample
        
        """

        # Change the dimensions from (b, c_in, z, y, x) to (b, c_in, x, y, z) so that the dimension order
        # is more intuitive
        U_x = U_x.permute(0, 1, 4, 3, 2)
        B_x = B_x.permute(0, 1, 4, 3, 2)

        x = torch.cat((U_x, B_x), dim=1)

        # ENCODER
        y_also_for_skip_connection_1 = self.encoder_convblock_1(x)
        y = self.downsample_1(y_also_for_skip_connection_1)

        y_also_for_skip_connection_2 = self.encoder_convblock_2(y)
        y = self.downsample_2(y_also_for_skip_connection_2)

        # BOTTLENECK
        y = self.bottleneck(y)

        # DECODER
        y = self.deconv_1(y)
        y = F.interpolate(y, size=y_also_for_skip_connection_2.shape[2:], mode="trilinear", align_corners=True)
        y = torch.cat((y, y_also_for_skip_connection_2), dim=1)
        y = self.decoder_convblock_1(y)

        y = self.deconv_2(y)
        y = F.interpolate(y, size=y_also_for_skip_connection_1.shape[2:], mode="trilinear", align_corners=True)
        y = torch.cat((y, y_also_for_skip_connection_1), dim=1)
        y = self.decoder_convblock_2(y)

        y = self.conv(y)

        # Change the dimensions from (b, c_out, x, y, z) to (b, c_out, z, y, x) so that the dimension order
        # is as needed by the C++ code that receives the output from forward()
        y = y.permute(0, 1, 4, 3, 2)

        return y

    
    def _block(in_channels, out_channels):
        return nn.Sequential(
            nn.Conv3d(in_channels, out_channels, kernel_size=3, padding=1),
            nn.ReLU(inplace=True),
            nn.Conv3d(out_channels, out_channels, kernel_size=3, padding=1),
            nn.ReLU(inplace=True)
            )
