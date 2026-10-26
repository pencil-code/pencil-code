import numpy as np
import pytest

from pencil.util import PencilArray

def test_pencilarray_noax():
	"""
	Test behaviour with no axes
	"""
	arr = np.arange(81).reshape((3,3,3,3))
	parr = PencilArray(arr)
	
	assert np.all(parr == arr)

@pytest.mark.xfail(reason="input validation not implemented?")
def test_pencilarray_axshape():
	"""
	Check that an error is raised if axis_order is specified with the wrong number of entries
	"""
	arr = np.arange(81).reshape((3,3,3,3))
	
	with pytest.raises(ValueError):
		_ = PencilArray(arr, axis_order=['x', 'y'])
	
	with pytest.raises(ValueError):
		_ = PencilArray(arr, axis_order=['x', 'y', 'z', 't', 'ivar'])

def test_pencilarray_axis_tracking():
	"""
	Is axis order correctly tracked through changes?
	"""
	arr = np.arange(24).reshape((2,3,4))
	parr = PencilArray(arr, axis_order=['z', 'y', 'x'])
	
	assert parr.axis_order == ('z', 'y', 'x')
	
	#indexing
	assert parr[:,1,:].axis_order == ('z', 'x')
	
	#slicing
	assert parr[:,:1,:].axis_order == ('z', 'y', 'x')
	
	#newaxis
	assert parr[:,None,:,:].axis_order == ('z', None, 'y', 'x')
	assert parr[:,np.newaxis,:,:].axis_order == ('z', None, 'y', 'x')


def test_pencilarray_axis_tracking_failure():
	"""
	Is axis_order set to None when it cannot be tracked?
	"""
	arr = np.arange(24).reshape((2,3,4))
	parr = PencilArray(arr, axis_order=['z', 'y', 'x'])
	
	#fancy indexing
	parr_2 = parr[[0,0,1]]
	assert parr_2.axis_order is None

@pytest.mark.xfail(reason="reshape not handled correctly")
def test_pencilarray_axis_tracking_failure_2():
	"""
	.reshape and transpose are not handled correctly; keep this test as xfail so that
	it is at least semi-documented that this doesn't work.
	"""
	arr = np.arange(24).reshape((2,3,4))
	parr = PencilArray(arr, axis_order=['z', 'y', 'x'])
	
	#reshape
	parr_3 = parr.reshape(6,4)
	with pytest.mark.xfail(reason="wrong implementation?"):
		assert parr_3.axis_order is None
	
	#transpose
	parr_4 = np.transpose(parr, axes=(1,0,2))
	assert parr_4.axis_order is None

def test_pencilarray_reorder():
	arr = np.arange(24).reshape((2,3,4))
	parr = PencilArray(arr, axis_order=['z', 'y', 'x'])
	
	parr_2 = parr.reorder('x', 'y', 'z')
	assert np.all(parr_2 == np.transpose(arr, axes=(2,1,0)))
