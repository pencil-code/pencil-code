# File containing functions for testing the kin-dynamo machine learning models.


import torch
import unittest

from kindynamo_lenet5forsegmentation import KinDynamo_LeNet5ForSegmentation
from kindynamo_tbnn import KinDynamo_TBNN
from kindynamo_unet import KinDynamo_UNet


class TestKinDynamoModels(unittest.TestCase):

    def setUp(self):

        torch.manual_seed(9)

        self.device = "cuda" if torch.cuda.is_available() else "cpu"

        self.batch_size = 1
        self.in_channels = 6
        self.x_size = 18
        self.y_size = 18
        self.z_size = 18

        self.out_channels = 3

        self.x = torch.rand((self.batch_size, self.in_channels, self.x_size, self.y_size, self.z_size), device=self.device)

        self.expected_output_size = torch.Size((self.batch_size, self.out_channels, self.x_size, self.y_size, self.z_size))


    ##### KinDynamo_LeNet5ForSegmentation #####

    def test_kindynamo_lenet5forsegmentation_output_size(self):
        """
        Test the forward() function output size of model KinDynamo_LeNet5ForSegmentation, defined in file kindynamo_lenet5forsegmentation.py.

        """

        # Test output size
        model = KinDynamo_LeNet5ForSegmentation().to(self.device)
        mean_residual_emf = model.forward(self.x)
        self.assertEqual(mean_residual_emf.size(), self.expected_output_size)



    def test_kindynamo_lenet5forsegmentation_output_values(self):
        """
        Test the forward() function output values of model KinDynamo_LeNet5ForSegmentation, defined in file kindynamo_lenet5forsegmentation.py.
        """

        # Test output values
        # TODO
        pass


    ##### KinDynamo_TBNN #####

    def test_kindynamo_tbnn_output_size(self):
        """
        Test the forward() function output size of model KinDynamo_TBNN, defined in file kindynamo_tbnn.py.

        """

        second_level_channels = 16
        spacing_val = 5

        # Test output size
        model = KinDynamo_TBNN(in_channels=self.in_channels, out_channels=self.out_channels, second_level_channels=second_level_channels, spacing_val=spacing_val).to(self.device)
        mean_residual_emf = model.forward(self.x)
        self.assertEqual(mean_residual_emf.size(), self.expected_output_size)


    def test_kindynamo_tbnn_output_values(self):
        """
        Test the forward() function output values of model KinDynamo_TBNN, defined in file kindynamo_tbnn.py.

        """

        # Test output values
        # TODO
        pass


    ##### KinDynamo_UNet #####

    def test_kindynamo_unet_output_size(self):
        """
        Test the forward() function output size of model KinDynamo_UNet, defined in file kindynamo_unet.py.

        """

        second_level_channels = 16

        # Test output size
        model = KinDynamo_UNet(in_channels=self.in_channels, out_channels=self.out_channels, second_level_channels=second_level_channels).to(self.device)
        mean_residual_emf = model.forward(self.x)
        self.assertEqual(mean_residual_emf.size(), self.expected_output_size)


    def test_kindynamo_output_values(self):
        """
        Test the forward() function output values of model KinDynamo_UNet, defined in file kindynamo_unet.py.

        """

        # Test output values
        # TODO
        pass


    if __name__ == "__main__":
        unittest.main()
