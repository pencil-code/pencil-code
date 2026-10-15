# File implementing a segmentation adaptation of the LeNet-5 for the kin-dynamo experiments.
# Concept for the original (classification) version from LeCun et al. (1998): "Gradient-Based Learning Applied to Document Recognition".


import torch
import torch.nn as nn
import torch.nn.functional as F

from collections import OrderedDict


# TODO: Create all missing class and method docstrings
# TODO: Complete all incomplete class and method docstrings (include both class or method descriptions and
#   parameter definitions)


class KinDynamo_LeNet5ForSegmentation(nn.Module):
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

   NOTE: LeNet5 was originally designed for image classification.
      However, because we wish to perform image segmentation, we remove
      LeNet5's fully connected layers and add deconvolutional layers.

   NOTE 2: The number of input channels is 3 here instead of 1, as LeNet-5 was
      designed for black-and-white images, but we are working with vector
      fields of 3-D magnetic fields and 3-D velocity fields

   NOTE 3: ReLU was not used in the original LeNet5 but is used here because
      it has proven useful in CNN architectures in the years since.

   NOTE 4: Max pooling is used instead of the original LeNet5's subsampling
      approach for simplicity, as it is available out-of-the-box from
      PyTorch. Average pooling could alternatively be explored.

   """

   def __init__(self):

      super(KinDynamo_LeNet5ForSegmentation, self).__init__()

      # First convolutional layer (and ReLU)
      self.conv_1 = nn.Conv3d(6, 6, kernel_size=3)
      self.relu_1 = nn.ReLU()

      # First downsampling layer
      self.downsample_1 = nn.MaxPool3d(2, stride=2)

      # Second convolutional layer (and ReLU)
      self.conv_2 = nn.Conv3d(6, 16, kernel_size=3)
      self.relu_2 = nn.ReLU()

      # Second downsampling layer
      self.downsample_2 = nn.MaxPool3d(2, stride=2)

      # Third convolutional layer (and ReLU)
      self.conv_3 = nn.Conv3d(16, 120, kernel_size=3)
      self.relu_3 = nn.ReLU()

      # ENCODER
      #####
      # DECODER

      # First deconvolutional layer (and ReLU)
      self.deconv_1 = nn.ConvTranspose3d(120, 16, kernel_size=3)
      self.relu_4 = nn.ReLU()

      # First upsampling layer
      # NOTE: Implemented in forward() using torch.nn.functional to support model inputs of varying size

      # Second deconvolutional layer (and ReLU)
      self.deconv_2 = nn.ConvTranspose3d(16, 6, kernel_size=3)
      self.relu_5 = nn.ReLU()

      # Second upsampling layer
      # NOTE: Implemented in forward() using torch.nn.functional to support model inputs of varying size

      # Third deconvolutional layer (NO ReLU)
      self.deconv_3 = nn.ConvTranspose3d(6, 3, kernel_size=3)


   def forward(self, U_x, B_x):
      """
      Perform a forward pass of the neural network. (I.e., perform inference on the input sample.)

      :param U_x: Mean velocity field values of the inputs sample
      
      :param B_x: Mean magnetic field values of the input sample

      """

      # Change the dimensions from (b, c_in, z, y, x) to (b, c_in, x, y, z) so that the dimension order
      # is more intuitive
      U_x = U_x.permute(0, 1, 4, 3, 2)
      B_x = B_x.permute(0, 1, 4, 3, 2)

      x = torch.cat((U_x, B_x), dim=1)

      y = self.conv_1(x)
      y = self.relu_1(y)

      pre_downsample_shape_1 = y.shape
      y = self.downsample_1(y)

      y = self.conv_2(y)
      y = self.relu_2(y)

      pre_downsample_shape_2 = y.shape
      y = self.downsample_2(y)

      y = self.conv_3(y)
      y = self.relu_3(y)

      # ENCODER
      #####
      # DECODER

      y = self.deconv_1(y)
      y = self.relu_4(y)

      y = F.interpolate(y, size=pre_downsample_shape_2[2:], mode="trilinear", align_corners=True)

      y = self.deconv_2(y)
      y = self.relu_5(y)

      y = F.interpolate(y, size=pre_downsample_shape_1[2:], mode="trilinear", align_corners=True)

      y = self.deconv_3(y)

      # Change the dimensions from (b, c_out, x, y, z) to (b, c_out, z, y, x) so that the dimension order
      # is as needed by the C++ code that receives the output from forward()
      y = y.permute(0, 1, 4, 3, 2)

      return y
