# **************************************************************************
# *
# * Authors:     J.M. de la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *              Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [2]
# *
# * [1] SciLifeLab, Stockholm University
# * [2] MRC Laboratory of Molecular Biology, MRC-LMB
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import pyworkflow.object as pwobj
from pyworkflow.utils import weakImport
from pwem.objects import EMObject, SetOfImages, Transform, CTFModel

with weakImport("tomo"):
    from tomo.objects import TiltImage, TomoAcquisition


class CryoDrgnParticles(EMObject):

    def __init__(self, filename=None, poses=None, ctfs=None,
                 dim=None, samplingRate=None, **kwargs):
        EMObject.__init__(self, **kwargs)

        self.filename = pwobj.String(filename)
        self.poses = pwobj.String(poses)
        self.ctfs = pwobj.String(ctfs)
        self.samplingRate = pwobj.Float(samplingRate)
        self.dim = pwobj.Integer(dim)
        self.ptcls = pwobj.Pointer()

    def __str__(self):
        return ('CryoDrgnParticles (%d x %d, %0.2f Å/px)'
                % (self.dim, self.dim, self.samplingRate))

    def getSamplingRate(self):
        return self.samplingRate.get()

    def getXDim(self):
        return self.dim.get()


class TiltParticle1(TiltImage):
    """ One tilted particle image from a tilt series. """

    def __init__(self, **kwargs):
        TiltImage.__init__(self, **kwargs)
        self._acquisition = TomoAcquisition()
        self._tsId = pwobj.String()
        self._volId = pwobj.Integer()
        self._coordinate = None  # 2D coord from a tilt series
        self._transform = None
        self._ctfModel = None

    def getTsId(self):
        return self._tsId.get()

    def setTsId(self, value):
        self._tsId.set(value)

    def getVolId(self):
        return self._volId.get()

    def setVolId(self, value):
        self._volId.set(value)

    def hasCoordinate(self):
        return self._coordinate is not None

    def setCoordinate(self, coordinate):
        self._coordinate = coordinate

    def getCoordinate(self):
        return self._coordinate

    def hasTransform(self):
        return self._transform is not None

    def getTransform(self) -> Transform:
        return self._transform

    def setTransform(self, newTransform):
        self._transform = newTransform

    def hasCTF(self):
        return self._ctfModel is not None

    def getCTF(self) -> CTFModel:
        """ Return the CTF model """
        return self._ctfModel

    def setCTF(self, newCTF):
        self._ctfModel = newCTF

    def __str__(self):
        dims = self.getDim()
        dimstr = dims or [None, None]
        return ('TiltParticle (%s x %s, %0.2f Å/px, '
                'tiltAngle: %0.1f, tsId: %s, volId: %d)'
                % (dimstr[0], dimstr[1], self._samplingRate.get(),
                   self.getTiltAngle(), self.getTsId(), self.getVolId()))


class SetOfTiltSeriesParticles(SetOfImages):
    """ Equivalent to a set of subtomograms, but in 2D. """
    ITEM_TYPE = TiltParticle1

    def __init__(self, **kwargs):
        SetOfImages.__init__(self, **kwargs)
        self.counter_ts = None
        self.counter_subtomo = None

    def getTSIds(self):
        """ Returns all TS_IDs involved in the set."""
        attr = "_tsId"
        tsIds = self.aggregate(["MAX"], attr, [attr])
        tsIds = [d[attr] for d in tsIds]
        return tsIds

    def getVolIds(self):
        """ Returns all volIds involved in the set."""
        attr = "_volId"
        volIds = self.aggregate(["MAX"], attr, [attr])
        volIds = [d[attr] for d in volIds]
        return volIds

    def __str__(self):
        if self.counter_ts is None:
            self.counter_ts = len(self.getTSIds())
        if self.counter_subtomo is None:
            self.counter_subtomo = len(self.getVolIds())

        return ('SetOfTiltSeriesParticles (%d items, %0.2f Å/px, '
                'from %d tilt series, %d subtomos)'
                % (self.getSize(), self._samplingRate.get(),
                   self.counter_ts, self.counter_subtomo))
