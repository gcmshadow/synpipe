#!/usr/bin/env python
# encoding: utf-8

from __future__ import print_function
from __future__ import absolute_import
from builtins import str
import numpy as np

import pyfits as fits

import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import lsst.afw.image as afwImage
import lsst.pex.config as afwConfig

from lsst.pipe.tasks.fakes import BaseFakeSourcesConfig, BaseFakeSourcesTask
from lsst.pex.exceptions import InvalidParameterError

from . import FakeSourceLib as fsl


class PositionStarFakesConfig(BaseFakeSourcesConfig):
    starList = afwConfig.Field(dtype=str,
                               doc="Catalog of stars with mags ra/dec")
    seed = afwConfig.Field(dtype=int, default=1,
                           doc="Seed for random number generator")


class PositionStarFakesTask(BaseFakeSourcesTask):
    ConfigClass = PositionStarFakesConfig

    def __init__(self, **kwargs):
        BaseFakeSourcesTask.__init__(self, **kwargs)
        print("RNG seed:", self.config.seed)
        self.rng = afwMath.Random(seed=self.config.seed)
        self.npRand = np.random.RandomState(self.config.seed)
        try:
            self.starData = fits.open(self.config.starList)[1].data
        except Exception:
            raise

    def run(self, exposure, background):

        self.log.info("Adding fake stars at real positions")
        psf = exposure.getPsf()
        psfBBox = psf.computeImage().getBBox()
        margin = max(psfBBox.getWidth(), psfBBox.getHeight())/2 + 1

        PARENT = afwImage.PARENT
        md = exposure.getMetadata()
        expBBox = exposure.getBBox(PARENT)
        wcs = exposure.getWcs()

        for istar, star in enumerate(self.starData):
            try:
                starident = star["ID"]
            except KeyError:
                starident = istar + 1

            try:
                flux = exposure.getCalib().getFlux(float(star['mag']))
            except KeyError:
                raise KeyError("No mag column in %s" % self.config.starList)

            try:
                starCoord = afwGeom.SpherePoint(star['RA'], star['DEC'], afwGeom.degrees)
            except KeyError:
                raise KeyError("No RA/DEC column in table".format(self.config.starList))

            starXY = wcs.skyToPixel(starCoord)
            bboxI = exposure.getBBox(PARENT)
            bboxI.grow(int(margin))
            if not bboxI.contains(afwGeom.Point2I(starXY)):
                continue

            try:
                starImage = psf.computeImage(starXY)
            except InvalidParameterError:
                # This means an image was computed in an area where there was no data
                # continue on to the next star. This is most likely to occur when inserting
                # into coadds
                logmsg = "Skipping fake {} because no input images present at point {}"
                self.log.info(logmsg.format(starident, starXY))
                continue

            starImage *= flux
            starBBox = starImage.getBBox(PARENT)

            # Check that we're within the larger exposure, otherwise crop
            if expBBox.contains(starBBox) is False:
                newBBox = starImage.getBBox(PARENT)
                newBBox.clip(expBBox)
                if newBBox.getArea() <= 0:
                    self.log.info("Skipping fake %d" % starident)
                    continue
                self.log.info("Cropping FAKE%d from %s to %s" % (starident,
                                                                 str(starBBox), str(newBBox)))
                starImage = starImage.Factory(starImage, newBBox, PARENT)
                starBBox = newBBox

            starMaskedImage = afwImage.MaskedImageF(starImage.convertF())

            starMaskedImage.getMask().set(self.bitmask)

            md.set("FAKE%s" % str(starident), "%.3f, %.3f" % (starXY.getX(),
                                                              starXY.getY()))
            self.log.info("Adding fake %s at: %.1f,%.1f" % (str(starident),
                                                            starXY.getX(),
                                                            starXY.getY()))

            maskedImage = exposure.getMaskedImage()
            BBox = starMaskedImage.getBBox(PARENT)
            subMaskedImage = maskedImage.Factory(exposure.getMaskedImage(),
                                                 BBox,
                                                 PARENT)
            subMaskedImage += starMaskedImage
