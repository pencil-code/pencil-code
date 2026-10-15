# File implementing a Tensor Basis Neural Network- (TBNN-) inspired neural network for the kin-dynamo experiments.
# Concept for the original TBNN from Ling et al. (2016): "Reynolds averaged turbulence modelling using deep neural networks with embedded invariance".
# The architecture for this TBNN-inspired neural network is also based on the U-Net architecture from Ronneberger et al. (2015): "U-Net: Convolutional Networks for Biomedical Image Segmentation".


import torch
import torch.nn as nn
import torch.nn.functional as F


# TODO: Install opt-einsum locally.
# TODO: Install opt-einsum on CSC server.
# TODO: Note use of opt-einsum package and reason for its use (for faster computation and/or less memory usage, per https://docs.pytorch.org/docs/2.12/generated/torch.einsum.html) in paper.

# TODO: Create all missing class and method docstrings
# TODO: Complete all incomplete class and method docstrings (include both class or method descriptions and
#   parameter definitions)

class KinDynamo_TBNN(nn.Module):
    """

    Physics-equation-based neural network for inferring mean residual EMF from mean magnetic field and
    mean velocity field.

    Equation (using Einstein notation) for mean residual EMF:

    Ebar_i = alpha_ij * Bbar_j + eta_ijk * dBbar_j/dx_k,

    where Ebar is the mean residual EMF, Bbar is the mean magnetic field, dBbar/dx is the partial
    derivative of the magnetic field with respect to x (the spatial x-dimension), and alpha and eta
    are the turbulent transport coefficients. Alpha and eta are, for our use case, which is the case of
    large-scale kinematic turbulent dynamos, dependent solely on the mean magnetic field and the mean
    mean velocity field.

    # TODO: Update the following description to reflect the updated architecture:
    The neural network is implemented such that the input follows two paths, one corresponding to the
    alpha_ij * Bbar_j part of the equation (herein dubbed the "alpha part") amd one corresponding to
    the eta_ijk * dBbar_j/dx_k part (herein dubbed the "eta part"), whose outputs are summed together
    elementwise to produce the final output of the network. This construction makes explicit the addition
    that takes place in the equation, meaning that the network does not need to learn to perform it, the
    intention being that reducing how much the network needs to learn will help it to learn the equation
    more effectively.

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

    """

    # def __init__(self, in_channels=6, out_channels=3, second_level_channels=16):
    def __init__(self, in_channels, out_channels, second_level_channels, spacing_val):

        super(KinDynamo_TBNN, self).__init__()

        self.in_channels = in_channels
        self.out_channels = out_channels
        self.spacing_val = spacing_val

        channels = second_level_channels

        # FOR ALPHA:

        # ENCODER
        self.alpha_encoder_convblock_1 = KinDynamo_TBNN._block(in_channels, channels)
        self.alpha_downsample_1 = nn.MaxPool3d(kernel_size=2, stride=2)

        self.alpha_encoder_convblock_2 = KinDynamo_TBNN._block(channels, channels * 2)
        self.alpha_downsample_2 = nn.MaxPool3d(kernel_size=2, stride=2)

        # BOTTLENECK
        self.alpha_bottleneck = KinDynamo_TBNN._block(channels * 2, channels * 4)

        # DECODER
        self.alpha_deconv_1 = nn.ConvTranspose3d(channels * 4, channels * 2, kernel_size=2, stride=2)
        # NOTE: Upsampling is implemented in forward() using torch.nn.functional to support model inputs
        #   of varying size.
        self.alpha_decoder_convblock_1 = KinDynamo_TBNN._block(channels * 4, channels * 2)

        self.alpha_deconv_2 = nn.ConvTranspose3d(channels * 2, channels, kernel_size=2, stride=2)
        # NOTE: Upsampling is implemented in forward() using torch.nn.functional to support model inputs
        #   of varying size.
        self.alpha_decoder_convblock_2 = KinDynamo_TBNN._block(channels * 2, channels)

        self.alpha_conv = nn.Conv3d(channels, out_channels*out_channels, kernel_size=1)
        

        # FOR ETA:

        self.eta_encoder_convblock_1 = KinDynamo_TBNN._block(in_channels, channels)
        self.eta_downsample_1 = nn.MaxPool3d(kernel_size=2, stride=2)

        self.eta_encoder_convblock_2 = KinDynamo_TBNN._block(channels, channels * 2)
        self.eta_downsample_2 = nn.MaxPool3d(kernel_size=2, stride=2)

        # BOTTLENECK
        self.eta_bottleneck = KinDynamo_TBNN._block(channels * 2, channels * 4)

        # DECODER
        self.eta_deconv_1 = nn.ConvTranspose3d(channels * 4, channels * 2, kernel_size=2, stride=2)
        # NOTE: Upsampling is implemented in forward() using torch.nn.functional to support model inputs
        #   of varying size.
        self.eta_decoder_convblock_1 = KinDynamo_TBNN._block(channels * 4, channels * 2)

        self.eta_deconv_2 = nn.ConvTranspose3d(channels * 2, channels, kernel_size=2, stride=2)
        # NOTE: Upsampling is implemented in forward() using torch.nn.functional to support model inputs
        #   of varying size.
        self.eta_decoder_convblock_2 = KinDynamo_TBNN._block(channels * 2, channels)

        self.eta_conv = nn.Conv3d(channels, out_channels*out_channels*out_channels, kernel_size=1) 


    def forward(self, U_x, B_x):
        """
        Perform a forward pass of the neural network. (I.e., perform inference on the input sample.)

        :param U_x: Mean velocity field values of the input sample
        
        :param B_x: Mean magnetic field values of the input sample

        :param spacing_val: Distance between adjacent grid points in the magnetic field and velocity field
            values. (Note that this is to be a singular value, so the spacing of the magnetic field grid
            points and of the velocity field grid points and in all dimensions should be the same.)
        
        """

        # Change the dimensions from (b, c_in, z, y, x) to (b, c_in, x, y, z) so that the dimension order
        # is more intuitive
        U_x = U_x.permute(0, 1, 4, 3, 2)
        B_x = B_x.permute(0, 1, 4, 3, 2)

        x = torch.cat((U_x, B_x), dim=1)

        # ALPHA

        # ENCODER
        alpha_y_also_for_skip_connection_1 = self.alpha_encoder_convblock_1(x)
        alpha_y = self.alpha_downsample_1(alpha_y_also_for_skip_connection_1)

        alpha_y_also_for_skip_connection_2 = self.alpha_encoder_convblock_2(alpha_y)
        alpha_y = self.alpha_downsample_2(alpha_y_also_for_skip_connection_2)

        # BOTTLENECK
        alpha_y = self.alpha_bottleneck(alpha_y)

        # DECODER
        alpha_y = self.alpha_deconv_1(alpha_y)
        alpha_y = F.interpolate(alpha_y, size=alpha_y_also_for_skip_connection_2.shape[2:], mode="trilinear", align_corners=True)
        alpha_y = torch.cat((alpha_y, alpha_y_also_for_skip_connection_2), dim=1)
        alpha_y = self.alpha_decoder_convblock_1(alpha_y)

        alpha_y = self.alpha_deconv_2(alpha_y)
        alpha_y = F.interpolate(alpha_y, size=alpha_y_also_for_skip_connection_1.shape[2:], mode="trilinear", align_corners=True)
        alpha_y = torch.cat((alpha_y, alpha_y_also_for_skip_connection_1), dim=1)
        alpha_y = self.alpha_decoder_convblock_2(alpha_y)

        alpha_y = self.alpha_conv(alpha_y)

        # Make alpha be a 2-D matrix at each grid point
        alpha_y = alpha_y.reshape(x.shape[0], self.out_channels, self.out_channels, x.shape[2], x.shape[3], x.shape[4])


        # ETA

        # ENCODER
        eta_y_also_for_skip_connection_1 = self.eta_encoder_convblock_1(x)
        eta_y = self.eta_downsample_1(eta_y_also_for_skip_connection_1)

        eta_y_also_for_skip_connection_2 = self.eta_encoder_convblock_2(eta_y)
        eta_y = self.eta_downsample_2(eta_y_also_for_skip_connection_2)

        # BOTTLENECK
        eta_y = self.eta_bottleneck(eta_y)

        # DECODER
        eta_y = self.eta_deconv_1(eta_y)
        eta_y = F.interpolate(eta_y, size=eta_y_also_for_skip_connection_2.shape[2:], mode="trilinear", align_corners=True)
        eta_y = torch.cat((eta_y, eta_y_also_for_skip_connection_2), dim=1)
        eta_y = self.eta_decoder_convblock_1(eta_y)

        eta_y = self.eta_deconv_2(eta_y)
        eta_y = F.interpolate(eta_y, size=eta_y_also_for_skip_connection_1.shape[2:], mode="trilinear", align_corners=True)
        eta_y = torch.cat((eta_y, eta_y_also_for_skip_connection_1), dim=1)
        eta_y = self.eta_decoder_convblock_2(eta_y)

        eta_y = self.eta_conv(eta_y)

        # Make eta be a 3-D matrix at each grid point
        eta_y = eta_y.reshape(x.shape[0], self.out_channels, self.out_channels, self.out_channels, x.shape[2], x.shape[3], x.shape[4])

        
        # MULTIPLY ALPHA AND MEAN MAGNETIC FIELD

        # Permute alpha to put the 2-D matrix at each grid point last
        alpha_y_permuted = torch.permute(alpha_y, (0, 3, 4, 5, 1, 2)) 

        # Permute B to put the vector at each grid point last
        B_x_permuted = torch.permute(B_x, (0, 2, 3, 4, 1))

        # Perform the multiplication
        alpha_times_B_y = torch.einsum("...ij,...j->...i", [alpha_y_permuted, B_x_permuted])

        # Permute product back to original dimension order
        alpha_times_B_y = torch.permute(alpha_times_B_y, (0, 4, 1, 2, 3))


        # MULTIPLY ETA AND DERIVATIVE OF MEAN MAGNETIC FIELD

        # Permute eta to put the 3-D matrix at each grid point last
        eta_y_permuted = torch.permute(eta_y, (0, 4, 5, 6, 1, 2, 3))

        # Calculate dB/dx and put it last
        dBdx = torch.stack(torch.gradient(B_x, dim=[-3, -2, -1], spacing=self.spacing_val), dim=-1)

        # Permute dB/dx to put the vector at each grid point second-to-last
        dBdx_permuted = torch.permute(dBdx, (0, 2, 3, 4, 1, 5))
        
        # Perform the multiplication
        eta_times_dBdx_y = torch.einsum("...ijk,...jk->...i", [eta_y_permuted, dBdx_permuted])

        # Permute product back to original dimension order
        eta_times_dBdx_y = torch.permute(eta_times_dBdx_y, (0, 4, 1, 2, 3))


        # ADD THE TWO MULTIPLIED PARTS

        y = torch.add(alpha_times_B_y, eta_times_dBdx_y)


        # Change the dimensions from (b, c_out, x, y, z) to (b, c_out, z, y, x) so that the dimension order
        # is as needed by the C++ code that receives the output from forward()
        y = y.permute(0, 1, 4, 3, 2)

        return y

    
    def get_spacing_val():
        return self.spacing_val


    def set_spacing_val(spacing_val):
        self.spacing_val = spacing_val


    def _block(in_channels, out_channels):
        return nn.Sequential(
            nn.Conv3d(in_channels, out_channels, kernel_size=3, padding=1),
            nn.ReLU(inplace=True),
            nn.Conv3d(out_channels, out_channels, kernel_size=3, padding=1),
            nn.ReLU(inplace=True)
            )
