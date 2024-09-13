#!/usr/bin/env python
# Authors: Chongzhi Zang
#
# This software is distributable under the terms of the GNU General
# Public License (GPL) v2, the text of which can be found at
# http://www.gnu.org/copyleft/gpl.html. Installing, importing or
# otherwise using this module constitutes acceptance of the terms of
# this License.
#
# Disclaimer
#
# This software is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# Comments and/or additions are welcome (send e-mail to:
# chongzhizang@gmail.com)

"""
This module contains classes of genome data, e.g. chromsomes
per species, the size of the chromosomes, etc.
"""

GenomeDataError = "Error in GenomeData class"

bg_number_chroms = 1
bg_length_of_chrom = 100000000

background_chroms = []
background_chrom_lengths = {}
for i in range(0, bg_number_chroms):
    background_chroms.append('chr' + str(i + 1))
    background_chrom_lengths['chr' + str(i + 1)] = bg_length_of_chrom

mm8_chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
              'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
              'chr18', 'chr19', 'chrX', 'chrY', 'chrM']

mm9_chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
              'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
              'chr18', 'chr19', 'chrX', 'chrY', 'chrM']

mm10_chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
               'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
               'chr18', 'chr19', 'chrX', 'chrY', 'chrM']

mm39_chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
               'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
               'chr18', 'chr19', 'chrX', 'chrY', 'chrM']

rn4_chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
              'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
              'chr18', 'chr19', 'chr20', 'chrX', 'chrM']

hg18_chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
               'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
               'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM']

hg19_chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
               'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
               'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM']

hg38_chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
               'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
               'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM']

sacCer1_chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
                  'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chrM']

sacCer3_chroms = {'chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrVI', 'chrVII', 'chrVIII', 'chrIX',
                         'chrX', 'chrXI', 'chrXII', 'chrXIII', 'chrXIV', 'chrXV', 'chrXVI'}

dm2_chroms = ['chr2h', 'chr2L', 'chr2R', 'chr3h', 'chr3L',
              'chr3R', 'chr4', 'chr4h', 'chrM', 'chrU', 'chrX', 'chrXh', 'chrYh']

dm3_chroms = ['chr2L', 'chr2LHet', 'chr2R', 'chr2RHet', 'chr3L', 'chr3LHet',
              'chr3R', 'chr3RHet', 'chr4', 'chrX', 'chrXHet', 'chrYHet', 'chrU', 'chrUextra', 'chrM']

dm6_chroms = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX', 'chrY']

pombe_chroms = ['chr1', 'chr2', 'chr3', 'mat']

tair8_chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5']

NC12_chroms = ['NC_026501.1', 'NC_026502.1', 'NC_026503.1', 'NC_026504.1', 'NC_026505.1', 'NC_026506.1', 'NC_026507.1',
               'NW_011929459.1', 'NW_011929460.1', 'NW_011929461.1', 'NW_011929462.1', 'NW_011929463.1', 'NW_011929464.1',
               'NW_011929465.1', 'NW_011929466.1', 'NW_011929467.1', 'NW_011929468.1', 'NW_011929469.1', 'NW_011929470.1',
               'NW_011929471.1', 'NC_026614.1']



mm8_chrom_lengths = {'chr1': 197069962, 'chr2': 181976762, 'chr3': 159872112,
                     'chr4': 155029701, 'chr5': 152003063, 'chr6': 149525685,
                     'chr7': 145134094, 'chr8': 132085098, 'chr9': 124000669,
                     'chr10': 129959148, 'chr11': 121798632, 'chr12': 120463159,
                     'chr13': 120614378, 'chr14': 123978870, 'chr15': 103492577,
                     'chr16': 98252459, 'chr17': 95177420, 'chr18': 90736837,
                     'chr19': 61321190, 'chrX': 165556469, 'chrY': 16029404,
                     'chrM': 16299}

mm9_chrom_lengths = {'chr1': 197195432, 'chr2': 181748087, 'chr3': 159599783,
                     'chr4': 155630120, 'chr5': 152537259, 'chr6': 149517037,
                     'chr7': 152524553, 'chr8': 131738871, 'chr9': 124076172,
                     'chr10': 129993255, 'chr11': 121843856, 'chr12': 121257530,
                     'chr13': 120284312, 'chr14': 125194864, 'chr15': 103494974,
                     'chr16': 98319150, 'chr17': 95272651, 'chr18': 90772031,
                     'chr19': 61342430, 'chrX': 166650296, 'chrY': 15902555,
                     'chrM': 16299}

mm10_chrom_lengths = {'chr1': 195471971, 'chr2': 182113224, 'chr3': 160039680,
                      'chr4': 156508116, 'chr5': 151834684, 'chr6': 149736546,
                      'chr7': 145441459, 'chr8': 129401213, 'chr9': 124595110,
                      'chr10': 130694993, 'chr11': 122082543, 'chr12': 120129022,
                      'chr13': 120421639, 'chr14': 124902244, 'chr15': 104043685,
                      'chr16': 98207768, 'chr17': 94987271, 'chr18': 90702639,
                      'chr19': 61431566, 'chrX': 171031299, 'chrY': 91744698,
                      'chrM': 16299}

mm39_chrom_lengths = {'chr1': 195154279, 'chr2': 181755017, 'chr3': 159745316,
                      'chr4': 156860686, 'chr5': 151758149, 'chr6': 149588044,
                      'chr7': 144995196, 'chr8': 130127694, 'chr9': 124359700,
                      'chr10': 130530862, 'chr11': 121973369, 'chr12': 120092757,
                      'chr13': 120883175, 'chr14': 125139656, 'chr15': 104073951,
                      'chr16': 98008968, 'chr17': 95294699, 'chr18': 90720763,
                      'chr19': 61420004, 'chrX': 169476592, 'chrY': 91455967,
                      'chrM': 16299}

rn4_chrom_lengths = {'chr1': 267910886, 'chr2': 258207540, 'chr3': 171063335,
                     'chr4': 187126005, 'chr5': 173096209, 'chr6': 147636619,
                     'chr7': 143002779, 'chr8': 129041809, 'chr9': 113440463,
                     'chr10': 110718848, 'chr11': 87759784, 'chr12': 46782294,
                     'chr13': 111154910, 'chr14': 112194335, 'chr15': 109758846,
                     'chr16': 90238779, 'chr17': 97296363, 'chr18': 87265094,
                     'chr19': 59218465, 'chr20': 55268282, 'chrX': 160699376,
                     'chrM': 16300}

hg18_chrom_lengths = {'chr1': 247249719, 'chr2': 242951149, 'chr3': 199501827,
                      'chr4': 191273063, 'chr5': 180857866, 'chr6': 170899992,
                      'chr7': 158821424, 'chr8': 146274826, 'chr9': 140273252,
                      'chr10': 135374737, 'chr11': 134452384, 'chr12': 132349534,
                      'chr13': 114142980, 'chr14': 106368585, 'chr15': 100338915,
                      'chr16': 88827254, 'chr17': 78774742, 'chr18': 76117153,
                      'chr19': 63811651, 'chr20': 62435964, 'chr21': 46944323,
                      'chr22': 49691432, 'chrX': 154913754, 'chrY': 57772954,
                      'chrM': 16571}

hg19_chrom_lengths = {'chr1': 249250621, 'chr2': 243199373, 'chr3': 198022430,
                      'chr4': 191154276, 'chr5': 180915260, 'chr6': 171115067,
                      'chr7': 159138663, 'chr8': 146364022, 'chr9': 141213431,
                      'chr10': 135534747, 'chr11': 135006516, 'chr12': 133851895,
                      'chr13': 115169878, 'chr14': 107349540, 'chr15': 102531392,
                      'chr16': 90354753, 'chr17': 81195210, 'chr18': 78077248,
                      'chr19': 59128983, 'chr20': 63025520, 'chr21': 48129895,
                      'chr22': 51304566, 'chrX': 155270560, 'chrY': 59373566,
                      'chrM': 16571}

hg38_chrom_lengths = {'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559,
                      'chr4': 190214555, 'chr5': 181538259, 'chr6': 170805979,
                      'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717,
                      'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
                      'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
                      'chr16': 90338345, 'chr17': 83257441, 'chr18': 80373285,
                      'chr19': 58617616, 'chr20': 64444167, 'chr21': 46709983,
                      'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415,
                      'chrM': 16569}

dm2_chrom_lengths = {'chr2h': 1694122, 'chr2L': 22407834, 'chr2R': 20766785,
                     'chr3h': 2955737, 'chr3L': 23771897, 'chr3R': 27905053,
                     'chr4h': 88110, 'chr4': 1281640, 'chrXh': 359526,
                     'chrX': 22224390, 'chrYh': 396896, 'chrM': 19517,
                     'chrU': 8724946}

dm3_chrom_lengths = {'chr2L': 23011544,
                     'chr2LHet': 368872,
                     'chr2R': 21146708,
                     'chr2RHet': 3288761,
                     'chr3L': 24543557,
                     'chr3LHet': 2555491,
                     'chr3R': 27905053,
                     'chr3RHet': 2517507,
                     'chr4': 1351857,
                     'chrX': 22422827,
                     'chrXHet': 204112,
                     'chrYHet': 347038,
                     'chrU': 10049037,
                     'chrUextra': 29004656,
                     'chrM': 19517}

dm6_chrom_lengths = {'chr2L': 23513712, 'chr2R': 25286936,
                     'chr3L': 28110227,'chr3R': 32079331,
                     'chr4': 1348131,
                     'chrX': 23542271, 'chrY': 3667352}

tair8_chrom_lengths = {'chr1': 30427671, 'chr2': 19698289, 'chr3': 23459830,
                       'chr4': 18585056, 'chr5': 26975502}

sacCer1_chrom_lengths = {'chr1': 230208, 'chr2': 813136, 'chr3': 316613,
                         'chr4': 1531914, 'chr5': 576869, 'chr6': 270148,
                         'chr7': 1090944, 'chr8': 562639, 'chr9': 439885,
                         'chr10': 745446, 'chr11': 666445, 'chr12': 1078173,
                         'chr13': 924430, 'chr14': 784328, 'chr15': 1091285,
                         'chr16': 948060, 'chrM': 85779}

sacCer3_chrom_lengths = {'chrI': 230218, 'chrII': 813184, 'chrIII': 316620,
                         'chrIV': 1531933, 'chrV': 576874, 'chrVI': 270161,
                         'chrVII': 1090940, 'chrVIII': 562643, 'chrIX': 439888,
                         'chrX': 745751, 'chrXI': 666816, 'chrXII': 1078177,
                         'chrXIII': 924431, 'chrXIV': 784333, 'chrXV': 1091291,
                         'chrXVI': 948066}

pombe_chrom_lengths = {'chr1': 5580032, 'chr2': 4541604, 'chr3': 2453783,
                       'mat': 41249}

NC12_chrom_lengths = {'NC_026501.1': 9798893, 'NC_026502.1': 4478683, 'NC_026503.1': 5274802,
                      'NC_026504.1': 6000761, 'NC_026505.1': 6436246, 'NC_026506.1': 4218384,
                      'NC_026507.1': 4255303, 'NW_011929459.1': 192308, 'NW_011929460.1': 142473,
                      'NW_011929461.1': 125404, 'NW_011929462.1': 31696, 'NW_011929463.1': 19714,
                      'NW_011929464.1': 13515, 'NW_011929465.1': 11565, 'NW_011929466.1': 9397,
                      'NW_011929467.1': 8983, 'NW_011929468.1': 6701, 'NW_011929469.1': 6309,
                      'NW_011929470.1': 4755, 'NW_011929471.1': 1646, 'NC_026614.1': 64840}



species_chroms = {'mm8': mm8_chroms,
                  'mm9': mm9_chroms,
                  'mm10': mm10_chroms,
                  'mm39': mm39_chroms,
                  'hg18': hg18_chroms,
                  'hg19': hg19_chroms,
                  'hg38': hg38_chroms,
                  'dm2': dm2_chroms,
                  'dm3': dm3_chroms,
                  'dm6': dm6_chroms,
                  'sacCer1': sacCer1_chroms,
                  'sacCer3': sacCer3_chroms,
                  'pombe': pombe_chroms,
                  'rn4': rn4_chroms,
                  'tair8': tair8_chroms,
                  'NC12': NC12_chroms,
                  'background': background_chroms}

species_chrom_lengths = {'mm8': mm8_chrom_lengths,
                         'mm9': mm9_chrom_lengths,
                         'mm10': mm10_chrom_lengths,
                         'mm39': mm39_chrom_lengths,
                         'hg18': hg18_chrom_lengths,
                         'hg19': hg19_chrom_lengths,
                         'hg38': hg38_chrom_lengths,
                         'dm2': dm2_chrom_lengths,
                         'dm3': dm3_chrom_lengths,
                         'dm6': dm6_chrom_lengths,
                         'sacCer1': sacCer1_chrom_lengths,
                         'sacCer3': sacCer3_chrom_lengths,
                         'pombe': pombe_chrom_lengths,
                         'rn4': rn4_chrom_lengths,
                         'tair8': tair8_chrom_lengths,
                         'NC12': NC12_chrom_lengths,
                         'background': background_chrom_lengths}
