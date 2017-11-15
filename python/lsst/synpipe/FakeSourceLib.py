#!/usr/bin/env python
# encoding: utf-8

import numpy as np

import lsst.afw.geom
import lsst.afw.math
import lsst.pex.config
import lsst.afw.image
import lsst.afw.cameraGeom
import lsst.pipe.base as pipeBase


"""
Helper functions for making fake sources
"""


class SkyMapIdContainer(pipeBase.DataIdContainer):
    """A version of lsst.pipe.base.DataIdContainer specialized for loading a skyMap
    in the make fake source catalog scripts. These scripts use the data id in a
    unique way, such that they only need a tract number. This class supports that
    use case and should not be used in any other contexts in the LSST stack.
    Required because butler.subset does not support only tract
    """

    def makeDataRefList(self, namespace):
        """Make self.refList from self.idList
        """

        for dataId in self.idList:
            if "tract" not in dataId:
                raise RuntimeError("id must specify which tract to process tract")
            # warn about unused options
            for key in dataId:
                if key != "tract":
                    namespace.log.warn("'{}' specified in --id is unused and will be ignored".format(key))
            addList = [dataId]

            self.refList += [namespace.butler.dataRef(datasetType="deepCoadd_skyMap", dataId=addId)
                             for addId in addList]

def cropFakeImage(fakeImage, expBBox):
    """
    Crops the Fake image to fit inside the exposure BBox
    Note that the bboxes need to have the correct offsets applied
    Args:
        fakeImage: fake image object
        expBBox:   bounding box for CCD exposure (integer type, BBoxI)
                   and with offsets applied

    Returns:
        New cropped fake image
    """
    fakeBBox = fakeImage.getBBox(lsst.afw.image.PARENT)

    if not expBBox.contains(fakeBBox):
        newBBox = fakeImage.getBBox(lsst.afw.image.PARENT)
        newBBox.clip(expBBox)
        fakeImage = fakeImage.Factory(fakeImage, newBBox,
                                      lsst.afw.image.PARENT)
        # TODO: finish this up


def addNoise(galImage, detector, rand_gen=None):
    """
    adds noise to the the image and returns a variance plane
    INPUT: image to add noise to
           detector where the image will be located, this sets the gain
    NOTE: this assumes float type images and will break if given doubles
    RETURN: a MaskedImageF with the image with additional noise and the
            variance plane
    giving the variance due to the object
    """
    # TODO: this is gaussian noise right now, probably good enough
    varImage = galImage.Factory(galImage, True)
    if rand_gen is None:
        rand_gen = np.random
    scale = np.sqrt(np.abs(varImage.getArray())) + 1e-12
    noiseArray = rand_gen.normal(loc=0.0,
                                 scale=scale,
                                 size=(galImage.getHeight(),
                                       galImage.getWidth()))
    noiseImage = lsst.afw.image.ImageF(noiseArray.astype(np.float32))
    galImage += noiseImage

    return lsst.afw.image.MaskedImageF(galImage, None, varImage)
