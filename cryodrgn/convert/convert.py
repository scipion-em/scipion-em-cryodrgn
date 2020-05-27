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
import os
import numpy as np
from collections import OrderedDict

import pyworkflow.utils as pwutils
import pwem
import pwem.convert.transformations as tfs
from pwem.emlib.image import ImageHandler

from .metadata import Table, Column


class WriterBase:
    """ Helper class to convert from Scipion SetOfImages subclasses
    into Relion>3.1 star files (and binaries if conversion needed).
    """
    def __init__(self, **kwargs):
        """
        Create a new instance with some configuration parameters.

        Keyword Args:
            rootDir: Specify a directory that will be used as "root"
                for setting the path values in the star file pointing
                to the binaries
            convertPolicy: By default, conversion of binary files will
                create symbolic links if no conversion is required.
                By passing convertPolicy=CONVERT_ALWAYS, it will force
                the conversion.
            useBaseName: (bool) By default the writer will use the id to
                generate shorter names. If this option is True, then
                the images base name will be used instead. This option
                might be useful in export protocols.

        """
        self._optics = None
        # Not used now
        #self.convertPolicy = kwargs.get('convertPolicy', self.CONVERT_IF_NEEDED)
        self.rootDir = kwargs.get('rootDir', None)
        self.outputDir = kwargs.get('outputDir', None)
        self.useBaseName = kwargs.get('useBaseName', False)
        self.extensions = kwargs.get('extensions', ['mrc'])
        self._ih = ImageHandler()  # used to convert images
        self._filesDict = {}  # used to map file names (converted or linked)
        self._dimensionality = 2
        self._imageSize = None

    def _createTableFromDict(self, rowDict):
        """ Helper function to create a Table instance from
        an input dict with keys as columns names and type
        the type of the values in the dict.
        """
        return Table(columns=[
            Column(k, type=type(v)) for k, v in rowDict.items()])

    def _ctfToRow(self, ctf, row):
        psd = ctf.getPsdFile()
        if psd:
            row['rlnCtfImage'] = psd
        dU, dV, dAngle = ctf.getDefocus()
        row['rlnDefocusU'] = dU
        row['rlnDefocusV'] = dV
        row['rlnCtfAstigmatism'] = abs(dU-dV)
        row['rlnDefocusAngle'] = dAngle
        row['rlnCtfFigureOfMerit'] = ctf.getFitQuality() or 0
        row['rlnCtfMaxResolution'] = ctf.getResolution() or 0

        phaseShift = ctf.getPhaseShift()

        if phaseShift is not None:
            row['rlnCtfPhaseShift'] = phaseShift


class Writer(WriterBase):
    """ Helper class to convert from Scipion SetOfImages subclasses
    into Relion>3.1 star files (and binaries if conversion needed).
    """
    def _getOpticsGroupNumber(self, img):
        """ Get the optics group number based on acquisition.
        Params:
            img: input image, movie, micrograph or particle
        """
        # Add now the new Optics Group stuff
        acq = img.getAcquisition()
        ogName = acq.opticsGroupName.get() or 'DefaultOpticsGroup'
        ps = img.getSamplingRate()

        if ogName not in self._optics:
            ogNumber = len(self._optics) + 1
            self._optics[ogName] = {
                'rlnOpticsGroupName': ogName,
                'rlnOpticsGroup': ogNumber,
                #'rlnMtfFileName': acq.mtfFile.get() or 'No-MTF',
                # FIXME: Check when we need to update the following
                'rlnMicrographOriginalPixelSize': ps,
                self._imgLabelPixelSize: ps,
                'rlnVoltage': acq.getVoltage(),
                'rlnSphericalAberration': acq.getSphericalAberration(),
                'rlnAmplitudeContrast': acq.getAmplitudeContrast(),
                'rlnBeamTiltX': acq.beamTiltX.get() or 0.,
                'rlnBeamTiltY': acq.beamTiltY.get() or 0.,
                'rlnImageDimensionality': self._dimensionality,
                'rlnImageSize': self._imageSize,
            }
            mtfFile = acq.mtfFile.get()
            if mtfFile is not None:
                self._optics[ogName]['rlnMtfFileName'] = mtfFile
        else:
            ogNumber = self._optics[ogName]['rlnOpticsGroup']

        return ogNumber

    def _setAttributes(self, obj, row, attributes):
        for attr in attributes:
            attrLabel = '_%s' % attributes
            if hasattr(obj, attrLabel):
                row[attr] = obj.getAttributeValue(attrLabel)

    def _align2DToRow(self, alignment, row):
        matrix = alignment.getMatrix()
        shifts = tfs.translation_from_matrix(matrix)
        shifts *= self._pixelSize
        angles = -np.rad2deg(tfs.euler_from_matrix(matrix, axes='szyz'))
        row['rlnOriginXAngst'], row['rlnOriginYAngst'] = shifts[:2]
        row['rlnAnglePsi'] = -(angles[0] + angles[2])

    def _alignProjToRow(self, alignment, row):
        matrix = np.linalg.inv(alignment.getMatrix())
        shifts = -tfs.translation_from_matrix(matrix)
        shifts *= self._pixelSize
        angles = -np.rad2deg(tfs.euler_from_matrix(matrix, axes='szyz'))
        row['rlnOriginXAngst'], row['rlnOriginYAngst'], row['rlnOriginZAngst'] = shifts
        row['rlnAngleRot'], row['rlnAngleTilt'], row['rlnAnglePsi'] = angles

    def _partToRow(self, part, row):
        row['rlnImageId'] = part.getObjId()

        # Add coordinate information
        coord = part.getCoordinate()
        if coord is not None:
            x, y = coord.getPosition()
            row['rlnCoordinateX'] = x
            row['rlnCoordinateY'] = y
            # Add some specify coordinate attributes
            self._setAttributes(coord, row, ['rlnClassNumber',
                                             'rlnAutopickFigureOfMerit',
                                             'rlnAnglePsi'])
            micName = coord.getMicName()
            if micName:
                row['rlnMicrographName'] = str(micName.replace(" ", ""))
            else:
                if coord.getMicId():
                    row['rlnMicrographName'] = str(coord.getMicId())

        index, fn = part.getLocation()
        if self.outputStack:
            row['rlnOriginalParticleName'] = locationToRelion(index, fn)
            index, fn = self._counter, self._relOutputStack
            if self._counter > 0:
                self._ih.convert(part, (index, self.outputStack))
        else:
            if self.outputDir is not None:
                fn = self._filesDict.get(fn, fn)

        row['rlnImageName'] = locationToRelion(index, fn)

        if self._setRandomSubset:
            row['rlnRandomSubset'] = part._rlnRandomSubset.get()

        # Set CTF values
        if self._setCtf:
            self._ctfToRow(part.getCTF(), row)

        # Set alignment if necessary
        if self._setAlign:
            self._setAlign(part.getTransform(), row)

        # Set additional labels if present
        self._setAttributes(part, row, self._extraLabels)

        # Add now the new Optics Group stuff
        row['rlnOpticsGroup'] = self._getOpticsGroupNumber(part)

        self._counter += 1

    def writeSetOfParticles(self, partsSet, starFile, **kwargs):
        # Process the first item and create the table based
        # on the generated columns
        self._imgLabelPixelSize = 'rlnImagePixelSize'

        self._optics = OrderedDict()
        partRow = OrderedDict()
        firstPart = partsSet.getFirstItem()

        # Convert binaries if required
        self.outputStack = kwargs.get('outputStack', None)
        if self.outputStack:
            self._relOutputStack = os.path.relpath(self.outputStack,
                                                   os.path.dirname(starFile))
        if self.outputDir is not None:
            forceConvert = kwargs.get('forceConvert', False)
            self._filesDict = convertBinaryFiles(partsSet, self.outputDir,
                                                 forceConvert=forceConvert)

        # Compute some flags from the first particle...
        # when flags are True, some operations will be applied to all particles
        self._preprocessImageRow = kwargs.get('preprocessImageRow', None)
        self._setRandomSubset = (kwargs.get('fillRandomSubset') and
                                 firstPart.hasAttribute('_rlnRandomSubset'))

        self._setCtf = kwargs.get('writeCtf', True) and firstPart.hasCTF()

        alignType = kwargs.get('alignType', partsSet.getAlignment())

        if alignType == pwem.ALIGN_2D:
            self._setAlign = self._align2DToRow
        elif alignType == pwem.ALIGN_PROJ:
            self._setAlign = self._alignProjToRow
        elif alignType == pwem.ALIGN_3D:
            raise Exception(
                "3D alignment conversion for Relion not implemented. "
                "It seems the particles were generated with an incorrect "
                "alignment type. You may either re-launch the protocol that "
                "generates the particles with angles or set 'Consider previous"
                " alignment?' to No")
        elif alignType == pwem.ALIGN_NONE:
            self._setAlign = None
        else:
            raise Exception("Invalid value for alignType: %s" % alignType)

        self._extraLabels = kwargs.get('extraLabels', [])
        self._extraLabels.extend(['rlnParticleSelectZScore',
                                  'rlnMovieFrameNumber'])
        self._postprocessImageRow = kwargs.get('postprocessImageRow', None)

        self._imageSize = firstPart.getXDim()
        self._pixelSize = firstPart.getSamplingRate() or 1.0

        self._counter = 0  # Mark first conversion as special one
        self._partToRow(firstPart, partRow)

        if self._postprocessImageRow:
            self._postprocessImageRow(firstPart, partRow)

        opticsTable = self._createTableFromDict(list(self._optics.values())[0])
        partsTable = self._createTableFromDict(partRow)
        partsTable.addRow(**partRow)

        with open(starFile, 'w') as f:
            # Write particles table
            f.write("# Star file generated with Scipion\n")
            f.write("# version 30001\n")
            # Write header first
            partsTable.writeStar(f, tableName='particles', writeRows=False)
            # Write all rows
            for part in partsSet:
                self._partToRow(part, partRow)
                if self._postprocessImageRow:
                    self._postprocessImageRow(part, partRow)
                partsTable.writeStarLine(f, partRow.values())

            # Write Optics at the end
            for opticsDict in self._optics.values():
                opticsTable.addRow(**opticsDict)
            f.write("\n# version 30001\n")
            opticsTable.writeStar(f, tableName='optics')


def locationToRelion(index, filename):
    """ Convert an index and filename location
    to a string with @ as expected in Relion.
    """
    if index != pwem.NO_INDEX:
        return "%06d@%s" % (index, filename)

    return filename


def convertBinaryFiles(imgSet, outputDir, extension='mrcs', forceConvert=False):
    """ Convert binary images files to a format read by Relion.
    Or create links if there is no need to convert the binary files.

    Params:
        imgSet: input image set to be converted.
        outputDir: where to put the converted file(s)
        extension: extension accepted by the program
        forceConvert: if True, the files will be converted and no root will be used
    Return:
        A dictionary with old-file as key and new-file as value
        If empty, not conversion was done.
    """
    filesDict = {}
    ih = ImageHandler()
    outputRoot = outputDir if forceConvert else os.path.join(outputDir, 'input')
    # Get the extension without the dot
    stackFiles = imgSet.getFiles()
    ext = pwutils.getExt(next(iter(stackFiles)))[1:]
    rootDir = pwutils.commonPath(list(stackFiles))

    def getUniqueFileName(fn, extension):
        """ Get an unique file for either link or convert files.
        It is possible that the base name overlap if they come
        from different runs. (like particles.mrcs after relion preprocess)
        """
        newFn = os.path.join(outputRoot, pwutils.replaceBaseExt(fn, extension))
        newRoot = pwutils.removeExt(newFn)

        values = filesDict.values()
        counter = 1

        while newFn in values:
            counter += 1
            newFn = '%s_%05d.%s' % (newRoot, counter, extension)

        return newFn

    def createBinaryLink(fn):
        """ Just create a link named .mrcs to Relion understand
        that it is a binary stack file and not a volume.
        """
        newFn = getUniqueFileName(fn, extension)
        if not os.path.exists(newFn):
            pwutils.createLink(fn, newFn)
            print("   %s -> %s" % (newFn, fn))
        return newFn

    def convertStack(fn):
        """ Convert from a format that is not read by Relion
        to an spider stack.
        """
        newFn = getUniqueFileName(fn, 'mrcs')
        ih.convertStack(fn, newFn)
        print("   %s -> %s" % (newFn, fn))
        return newFn

    def replaceRoot(fn):
        """ Link create to the root folder, so just replace that
        in the name, no need to do anything else.
        """
        return fn.replace(rootDir, outputRoot)

    if forceConvert:
        mapFunc = convertStack
    elif ext == extension:
        print("convertBinaryFiles: creating soft links.")
        print("   Root: %s -> %s" % (outputRoot, rootDir))
        mapFunc = replaceRoot
        # FIXME: There is a bug in pwutils.createLink when input is a single folder
        # pwutils.createLink(rootDir, outputRoot)
        # relativeOutput = os.path.join(os.path.relpath(rootDir, outputRoot), rootDir)
        # If the rootDir is a prefix in the outputRoot (usually Runs)
        # we need to prepend that basename to make the link works
        if rootDir in outputRoot:
            relativeOutput = os.path.join(os.path.relpath(rootDir, outputRoot),
                                          os.path.basename(rootDir))
        else:
            relativeOutput = os.path.relpath(rootDir,
                                             os.path.dirname(outputRoot))
        if not os.path.exists(outputRoot):
            os.symlink(relativeOutput, outputRoot)
    elif ext == 'mrc' and extension == 'mrcs':
        print("convertBinaryFiles: creating soft links (mrcs -> mrc).")
        mapFunc = createBinaryLink
    elif ext.endswith('hdf'):  # assume eman .hdf format
        print("convertBinaryFiles: converting stacks. (%s -> %s)"
              % (extension, ext))
        mapFunc = convertStack
    else:
        mapFunc = None

    if mapFunc is not None:
        pwutils.makePath(outputRoot)
        for fn in stackFiles:
            newFn = mapFunc(fn)  # convert or link
            filesDict[fn] = newFn  # map new filename

    return filesDict
