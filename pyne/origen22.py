from __future__ import print_function, division
try:
    from future_builtins import map, zip
except ImportError:
    pass
import os
import re
import sys
from io import StringIO
from collections import Mapping
from copy import deepcopy
from itertools import chain
from warnings import warn
from pyne.utils import QAWarning

import numpy as np

from pyne import data
from pyne import nucname
from pyne.xs import cache
from pyne import decay_tape9
from pyne.material import Material, from_atom_frac

if sys.version_info[0] > 2:
    basestring = str
    unicode = str

warn(__name__ + " is not yet QA compliant.", QAWarning)

BASE_TAPE9 = os.path.join(os.path.dirname(__file__), 'base_tape9.inp')

ACTIVATION_PRODUCT_NUCS = frozenset([10010000,
    10020000,  10030000,  10040000,  20030000,  20040000,  20060000,  30060000,
    30070000,  30080000,  40080000,  40090000,  40100000,  40110000,  50100000,
    50110000,  50120000,  60120000,  60130000,  60140000,  60150000,  70130000,
    70140000,  70150000,  70160000,  80160000,  80170000,  80180000,  80190000,
    90190000,  90200000,  100200000, 100210000, 100220000, 100230000, 110220000,
    110230000, 110240000, 110240001, 110250000, 120240000, 120250000, 120260000,
    120270000, 120280000, 130270000, 130280000, 130290000, 130300000, 140280000,
    140290000, 140300000, 140310000, 140320000, 150310000, 150320000, 150330000,
    150340000, 160320000, 160330000, 160340000, 160350000, 160360000, 160370000,
    170350000, 170360000, 170370000, 170380000, 170380001, 180360000, 180370000,
    180380000, 180390000, 180400000, 180410000, 180420000, 190390000, 190400000,
    190410000, 190420000, 190430000, 190440000, 200400000, 200410000, 200420000,
    200430000, 200440000, 200450000, 200460000, 200470000, 200480000, 200490000,
    210450000, 210460000, 210460001, 210470000, 210480000, 210490000, 210500000,
    220460000, 220470000, 220480000, 220490000, 220500000, 220510000, 230490000,
    230500000, 230510000, 230520000, 230530000, 230540000, 240500000, 240510000,
    240520000, 240530000, 240540000, 240550000, 250540000, 250550000, 250560000,
    250570000, 250580000, 260540000, 260550000, 260560000, 260570000, 260580000,
    260590000, 270580000, 270580001, 270590000, 270600000, 270600001, 270610000,
    270620000, 280580000, 280590000, 280600000, 280610000, 280620000, 280630000,
    280640000, 280650000, 280660000, 290620000, 290630000, 290640000, 290650000,
    290660000, 290670000, 300630000, 300640000, 300650000, 300660000, 300670000,
    300680000, 300690000, 300690001, 300700000, 300710000, 300710001, 310690000,
    310700000, 310710000, 310720000, 310720001, 320700000, 320710000, 320710001,
    320720000, 320730000, 320740000, 320750000, 320750001, 320760000, 320770000,
    320770001, 330750000, 330760000, 330770000, 340740000, 340750000, 340760000,
    340770000, 340770001, 340780000, 340790000, 340790001, 340800000, 340810000,
    340810001, 340820000, 340830000, 340830001, 350790000, 350800000, 350800001,
    350810000, 350820000, 350820001, 350830000, 360780000, 360790000, 360790001,
    360800000, 360810000, 360810001, 360820000, 360830000, 360830001, 360840000,
    360850000, 360850001, 360860000, 360870000, 360880000, 370850000, 370860000,
    370860001, 370870000, 370880000, 370890000, 380840000, 380850000, 380850001,
    380860000, 380870000, 380870001, 380880000, 380890000, 380900000, 380910000,
    380930000, 390890000, 390890001, 390900000, 390900001, 390910000, 390920000,
    390930000, 390940000, 390960000, 400890000, 400900000, 400910000, 400920000,
    400930000, 400940000, 400950000, 400960000, 400970000, 410910000, 410920000,
    410930000, 410930001, 410940000, 410950000, 410950001, 410960000, 410970000,
    410970001, 410980000, 411000000, 420920000, 420930000, 420930001, 420940000,
    420950000, 420960000, 420970000, 420980000, 420990000, 421000000, 421010000,
    430970000, 430970001, 430980000, 430990000, 431000000, 431010000, 440960000,
    440970000, 440980000, 440990000, 441000000, 441010000, 441020000, 441030000,
    441040000, 441050000, 441060000, 441070000, 451020000, 451030000, 451040000,
    451040001, 451050000, 451050001, 451060000, 451060001, 451070000, 461020000,
    461030000, 461040000, 461050000, 461060000, 461070000, 461070001, 461080000,
    461090000, 461090001, 461100000, 461110000, 461110001, 471060000, 471070000,
    471080000, 471080001, 471090000, 471090001, 471100000, 471100001, 471110000,
    471110001, 471120000, 481060000, 481070000, 481080000, 481090000, 481100000,
    481110000, 481110001, 481120000, 481130000, 481140000, 481150000, 481150001,
    481160000, 481170000, 481170001, 481190000, 481210000, 491130000, 491130001,
    491140000, 491140001, 491150000, 491160000, 491160001, 491170000, 491170001,
    491180000, 491190000, 491190001, 491200000, 491200001, 491210000, 501120000,
    501130000, 501130001, 501140000, 501150000, 501160000, 501170000, 501170001,
    501180000, 501190000, 501190001, 501200000, 501210000, 501210001, 501220000,
    501230000, 501230001, 501240000, 501250000, 501250001, 511210000, 511220000,
    511220001, 511230000, 511240000, 511240001, 511250000, 511260000, 511260001,
    521200000, 521210000, 521210001, 521220000, 521230000, 521230001, 521240000,
    521250000, 521250001, 521260000, 521270000, 521270001, 521280000, 521290000,
    521290001, 521300000, 521310000, 521310001, 531250000, 531260000, 531270000,
    531280000, 531290000, 531300000, 531300001, 531310000, 531320000, 541240000,
    541250000, 541250001, 541260000, 541270000, 541270001, 541280000, 541290000,
    541290001, 541300000, 541310000, 541310001, 541320000, 541330000, 541330001,
    541340000, 541350000, 541350001, 541360000, 541370000, 551310000, 551320000,
    551330000, 551340000, 551340001, 551350000, 551360000, 551370000, 551380000,
    561300000, 561310000, 561310001, 561320000, 561330000, 561330001, 561340000,
    561350000, 561350001, 561360000, 561360001, 561370000, 561370001, 561380000,
    561390000, 561400000, 561410000, 571370000, 571380000, 571390000, 571400000,
    571410000, 581360000, 581370000, 581370001, 581380000, 581390000, 581390001,
    581400000, 581410000, 581420000, 581430000, 581440000, 581450000, 591410000,
    591420000, 591420001, 591430000, 591440000, 591450000, 601420000, 601430000,
    601440000, 601450000, 601460000, 601470000, 601480000, 601490000, 601500000,
    601510000, 611450000, 611470000, 611480000, 611480001, 611490000, 611500000,
    611510000, 611520000, 621440000, 621450000, 621460000, 621470000, 621480000,
    621490000, 621500000, 621510000, 621520000, 621530000, 621540000, 621550000,
    631510000, 631520000, 631520001, 631530000, 631540000, 631550000, 631560000,
    641520000, 641530000, 641540000, 641550000, 641550001, 641560000, 641570000,
    641580000, 641590000, 641600000, 641610000, 641620000, 651570000, 651590000,
    651600000, 651610000, 651620000, 661560000, 661570000, 661580000, 661590000,
    661600000, 661610000, 661620000, 661630000, 661640000, 661650000, 661650001,
    661660000, 671630000, 671650000, 671660000, 671660001, 681620000, 681630000,
    681640000, 681650000, 681660000, 681670000, 681670001, 681680000, 681690000,
    681700000, 681710000, 681720000, 691690000, 691700000, 691700001, 691710000,
    691720000, 691730000, 701680000, 701690000, 701700000, 701710000, 701720000,
    701730000, 701740000, 701750000, 701750001, 701760000, 701770000, 711750000,
    711760000, 711760001, 711770000, 711770001, 721740000, 721750000, 721760000,
    721770000, 721780000, 721780001, 721790000, 721790001, 721800000, 721800001,
    721810000, 721820000, 731800000, 731810000, 731820000, 731820001, 731830000,
    741800000, 741810000, 741820000, 741830000, 741830001, 741840000, 741850000,
    741850001, 741860000, 741870000, 741880000, 741890000, 751850000, 751860000,
    751870000, 751880000, 751880001, 751890000, 761840000, 761850000, 761860000,
    761870000, 761880000, 761890000, 761900000, 761900001, 761910000, 761910001,
    761920000, 761930000, 761940000, 771910000, 771920000, 771920001, 771930000,
    771940000, 771940001, 781900000, 781910000, 781920000, 781930000, 781930001,
    781940000, 781950000, 781950001, 781960000, 781970000, 781970001, 781980000,
    781990000, 781990001, 791970000, 791980000, 791990000, 792000000, 801960000,
    801970000, 801970001, 801980000, 801990000, 801990001, 802000000, 802010000,
    802020000, 802030000, 802040000, 802050000, 812030000, 812040000, 812050000,
    812060000, 822040000, 822050000, 822060000, 822070000, 822080000, 822090000,
    832080000, 832090000, 832100000, 832100001, 832110000, 842100000, 842110000,
    842110001])
"""Set of activation product nuclides in id form."""

ACTINIDE_AND_DAUGHTER_NUCS = frozenset([
    20040000,  812060000, 812070000, 812080000, 812090000, 822060000,
    822070000, 822080000, 822090000, 822100000, 822110000, 822120000, 822140000,
    832080000, 832090000, 832100000, 832100001, 832110000, 832120000, 832130000,
    832140000, 842100000, 842110000, 842110001, 842120000, 842130000, 842140000,
    842150000, 842160000, 842180000, 852170000, 862180000, 862190000, 862200000,
    862220000, 872210000, 872230000, 882220000, 882230000, 882240000, 882250000,
    882260000, 882280000, 892250000, 892270000, 892280000, 902260000, 902270000,
    902280000, 902290000, 902300000, 902310000, 902320000, 902330000, 902340000,
    912310000, 912320000, 912330000, 912340000, 912340001, 912350000, 922300000,
    922310000, 922320000, 922330000, 922340000, 922350000, 922360000, 922370000,
    922380000, 922390000, 922400000, 922410000, 932350000, 932360000, 932360001,
    932370000, 932380000, 932390000, 932400000, 932400001, 932410000, 942360000,
    942370000, 942380000, 942390000, 942400000, 942410000, 942420000, 942430000,
    942440000, 942450000, 942460000, 952390000, 952400000, 952410000, 952420000,
    952420001, 952430000, 952440000, 952440001, 952450000, 952460000, 962410000,
    962420000, 962430000, 962440000, 962450000, 962460000, 962470000, 962480000,
    962490000, 962500000, 962510000, 972490000, 972500000, 972510000, 982490000,
    982500000, 982510000, 982520000, 982530000, 982540000, 982550000, 992530000,
    992540000, 992540001, 992550000])
"""Set of actinide & daughter nuclides in id form."""

FISSION_PRODUCT_NUCS = frozenset([
    10030000,  30060000,  30070000,  40090000,  40100000,  60140000,  270720000,
    270730000, 270740000, 270750000, 280660000, 280720000, 280730000, 280740000,
    280750000, 280760000, 280770000, 280780000, 290660000, 290670000, 290720000,
    290730000, 290740000, 290750000, 290760000, 290770000, 290780000, 290790000,
    290800000, 290810000, 300660000, 300670000, 300680000, 300690000, 300690001,
    300700000, 300710000, 300710001, 300720000, 300730000, 300740000, 300750000,
    300760000, 300770000, 300780000, 300790000, 300800000, 300810000, 300820000,
    300830000, 310690000, 310700000, 310710000, 310720000, 310730000, 310740000,
    310750000, 310760000, 310770000, 310780000, 310790000, 310800000, 310810000,
    310820000, 310830000, 310840000, 310850000, 320700000, 320710000, 320710001,
    320720000, 320730000, 320730001, 320740000, 320750000, 320750001, 320760000,
    320770000, 320770001, 320780000, 320790000, 320800000, 320810000, 320820000,
    320830000, 320840000, 320850000, 320860000, 320870000, 320880000, 330750000,
    330760000, 330770000, 330780000, 330790000, 330800000, 330810000, 330820000,
    330820001, 330830000, 330840000, 330850000, 330860000, 330870000, 330880000,
    330890000, 330900000, 340760000, 340770000, 340770001, 340780000, 340790000,
    340790001, 340800000, 340810000, 340810001, 340820000, 340830000, 340830001,
    340840000, 340850000, 340850001, 340860000, 340870000, 340880000, 340890000,
    340900000, 340910000, 340920000, 340930000, 350790000, 350790001, 350800000,
    350800001, 350810000, 350820000, 350820001, 350830000, 350840000, 350840001,
    350850000, 350860000, 350860001, 350870000, 350880000, 350890000, 350900000,
    350910000, 350920000, 350930000, 350940000, 350950000, 350960000, 360790000,
    360800000, 360810000, 360810001, 360820000, 360830000, 360830001, 360840000,
    360850000, 360850001, 360860000, 360870000, 360880000, 360890000, 360900000,
    360910000, 360920000, 360930000, 360940000, 360950000, 360960000, 360970000,
    360980000, 370850000, 370860000, 370860001, 370870000, 370880000, 370890000,
    370900000, 370900001, 370910000, 370920000, 370930000, 370940000, 370950000,
    370960000, 370970000, 370980000, 370990000, 371000000, 371010000, 380860000,
    380870000, 380870001, 380880000, 380890000, 380900000, 380910000, 380920000,
    380930000, 380940000, 380950000, 380960000, 380970000, 380980000, 380990000,
    381000000, 381010000, 381020000, 381030000, 381040000, 390890000, 390890001,
    390900000, 390900001, 390910000, 390910001, 390920000, 390930000, 390940000,
    390950000, 390960000, 390970000, 390980000, 390990000, 391000000, 391010000,
    391020000, 391030000, 391040000, 391050000, 391060000, 391070000, 400900000,
    400900001, 400910000, 400920000, 400930000, 400940000, 400950000, 400960000,
    400970000, 400980000, 400990000, 401000000, 401010000, 401020000, 401030000,
    401040000, 401050000, 401060000, 401070000, 401080000, 401090000, 410910000,
    410920000, 410930000, 410930001, 410940000, 410940001, 410950000, 410950001,
    410960000, 410970000, 410970001, 410980000, 410980001, 410990000, 410990001,
    411000000, 411000001, 411010000, 411020000, 411030000, 411040000, 411050000,
    411060000, 411070000, 411080000, 411090000, 411100000, 411110000, 411120000,
    420950000, 420960000, 420970000, 420980000, 420990000, 421000000, 421010000,
    421020000, 421030000, 421040000, 421050000, 421060000, 421070000, 421080000,
    421090000, 421100000, 421110000, 421120000, 421130000, 421140000, 421150000,
    430980000, 430990000, 430990001, 431000000, 431010000, 431020000, 431020001,
    431030000, 431040000, 431050000, 431060000, 431070000, 431080000, 431090000,
    431100000, 431110000, 431120000, 431130000, 431140000, 431150000, 431160000,
    431170000, 431180000, 440990000, 441000000, 441010000, 441020000, 441030000,
    441040000, 441050000, 441060000, 441070000, 441080000, 441090000, 441100000,
    441110000, 441120000, 441130000, 441140000, 441150000, 441160000, 441170000,
    441180000, 441190000, 441200000, 451020000, 451030000, 451030001, 451040000,
    451040001, 451050000, 451050001, 451060000, 451060001, 451070000, 451080000,
    451080001, 451090000, 451090001, 451100000, 451100001, 451110000, 451120000,
    451130000, 451140000, 451150000, 451160000, 451170000, 451180000, 451190000,
    451200000, 451210000, 451220000, 451230000, 461020000, 461040000, 461050000,
    461060000, 461070000, 461070001, 461080000, 461090000, 461090001, 461100000,
    461110000, 461110001, 461120000, 461130000, 461140000, 461150000, 461160000,
    461170000, 461180000, 461190000, 461200000, 461210000, 461220000, 461230000,
    461240000, 461250000, 461260000, 471060000, 471070000, 471080000, 471080001,
    471090000, 471090001, 471100000, 471100001, 471110000, 471110001, 471120000,
    471130000, 471130001, 471140000, 471150000, 471150001, 471160000, 471160001,
    471170000, 471170001, 471180000, 471180001, 471190000, 471200000, 471210000,
    471220000, 471230000, 471240000, 471250000, 471260000, 471270000, 471280000,
    481080000, 481090000, 481100000, 481110000, 481110001, 481120000, 481130000,
    481130001, 481140000, 481150000, 481150001, 481160000, 481170000, 481170001,
    481180000, 481190000, 481190001, 481200000, 481210000, 481220000, 481230000,
    481240000, 481250000, 481260000, 481270000, 481280000, 481290000, 481300000,
    481310000, 481320000, 491130000, 491130001, 491140000, 491140001, 491150000,
    491150001, 491160000, 491160001, 491170000, 491170001, 491180000, 491180001,
    491190000, 491190001, 491200000, 491200001, 491210000, 491210001, 491220000,
    491220001, 491230000, 491230001, 491240000, 491250000, 491250001, 491260000,
    491270000, 491270001, 491280000, 491290000, 491300000, 491310000, 491320000,
    491330000, 491340000, 501140000, 501150000, 501160000, 501170000, 501170001,
    501180000, 501190000, 501190001, 501200000, 501210000, 501210001, 501220000,
    501230000, 501230001, 501240000, 501250000, 501250001, 501260000, 501270000,
    501270001, 501280000, 501290000, 501290001, 501300000, 501310000, 501320000,
    501330000, 501340000, 501350000, 501360000, 511210000, 511220000, 511220001,
    511230000, 511240000, 511240001, 511250000, 511260000, 511260001, 511270000,
    511280000, 511280001, 511290000, 511300000, 511300001, 511310000, 511320000,
    511320001, 511330000, 511340000, 511340001, 511350000, 511360000, 511370000,
    511380000, 511390000, 521220000, 521230000, 521230001, 521240000, 521250000,
    521250001, 521260000, 521270000, 521270001, 521280000, 521290000, 521290001,
    521300000, 521310000, 521310001, 521320000, 521330000, 521330001, 521340000,
    521350000, 521360000, 521370000, 521380000, 521390000, 521400000, 521410000,
    521420000, 531270000, 531280000, 531290000, 531300000, 531300001, 531310000,
    531320000, 531330000, 531330001, 531340000, 531340001, 531350000, 531360000,
    531360001, 531370000, 531380000, 531390000, 531400000, 531410000, 531420000,
    531430000, 531440000, 531450000, 541260000, 541270000, 541280000, 541290000,
    541290001, 541300000, 541310000, 541310001, 541320000, 541330000, 541330001,
    541340000, 541340001, 541350000, 541350001, 541360000, 541370000, 541380000,
    541390000, 541400000, 541410000, 541420000, 541430000, 541440000, 541450000,
    541460000, 541470000, 551320000, 551330000, 551340000, 551340001, 551350000,
    551350001, 551360000, 551370000, 551380000, 551380001, 551390000, 551400000,
    551410000, 551420000, 551430000, 551440000, 551450000, 551460000, 551470000,
    551480000, 551490000, 551500000, 561320000, 561330000, 561340000, 561350000,
    561350001, 561360000, 561360001, 561370000, 561370001, 561380000, 561390000,
    561400000, 561410000, 561420000, 561430000, 561440000, 561450000, 561460000,
    561470000, 561480000, 561490000, 561500000, 561510000, 561520000, 571380000,
    571390000, 571400000, 571410000, 571420000, 571430000, 571440000, 571450000,
    571460000, 571470000, 571480000, 571490000, 571500000, 571510000, 571520000,
    571530000, 571540000, 571550000, 581390000, 581400000, 581410000, 581420000,
    581430000, 581440000, 581450000, 581460000, 581470000, 581480000, 581490000,
    581500000, 581510000, 581520000, 581530000, 581540000, 581550000, 581560000,
    581570000, 591390000, 591400000, 591410000, 591420000, 591420001, 591430000,
    591440000, 591440001, 591450000, 591460000, 591470000, 591480000, 591490000,
    591500000, 591510000, 591520000, 591530000, 591540000, 591550000, 591560000,
    591570000, 591580000, 591590000, 601410000, 601420000, 601430000, 601440000,
    601450000, 601460000, 601470000, 601480000, 601490000, 601500000, 601510000,
    601520000, 601530000, 601540000, 601550000, 601560000, 601570000, 601580000,
    601590000, 601600000, 601610000, 611450000, 611460000, 611470000, 611480000,
    611480001, 611490000, 611500000, 611510000, 611520000, 611520001, 611530000,
    611540000, 611540001, 611550000, 611560000, 611570000, 611580000, 611590000,
    611600000, 611610000, 611620000, 621450000, 621460000, 621470000, 621480000,
    621490000, 621500000, 621510000, 621520000, 621530000, 621540000, 621550000,
    621560000, 621570000, 621580000, 621590000, 621600000, 621610000, 621620000,
    621630000, 621640000, 621650000, 631490000, 631500000, 631510000, 631520000,
    631520001, 631530000, 631540000, 631550000, 631560000, 631570000, 631580000,
    631590000, 631600000, 631610000, 631620000, 631630000, 631640000, 631650000,
    641520000, 641530000, 641540000, 641550000, 641550001, 641560000, 641570000,
    641580000, 641590000, 641600000, 641610000, 641620000, 641630000, 641640000,
    641650000, 651590000, 651600000, 651610000, 651620000, 651620001, 651630000,
    651630001, 651640000, 651650000, 661600000, 661610000, 661620000, 661630000,
    661640000, 661650000, 661650001, 661660000, 671650000, 671660000, 671660001,
    681660000, 681670000, 681670001, 681680000, 681690000, 681700000, 681710000,
    681720000, 691690000, 691700000, 691700001, 691710000, 691720000, 701680000,
    701690000, 701700000, 701710000, 701720000])
"""Set of fission product nuclides in id form."""

NUCS = ACTIVATION_PRODUCT_NUCS | ACTINIDE_AND_DAUGHTER_NUCS | FISSION_PRODUCT_NUCS
"""Set of all known nuclides."""

DECAY_FIELDS = ('half_life', 'frac_beta_minus_x',
    'frac_beta_plus_or_electron_capture',
    'frac_beta_plus_or_electron_capture_x',
    'frac_alpha', 'frac_isomeric_transition', 'frac_spont_fiss', 'frac_beta_n',
    'recoverable_energy', 'frac_natural_abund', 'inhilation_concentration',
    'ingestion_concentration')
"""The decay data keys in a tape9 dictionary."""

XSFPY_FIELDS = ('sigma_gamma', 'sigma_2n', 'sigma_gamma_x', 'sigma_2n_x',
    'fiss_yields_present')
"""The cross section and fission product yield data keys in a tape9 dictionary.
"""

ACTIVATION_PRODUCT_FIELDS = ('sigma_3n', 'sigma_p')
"""The cross section data keys for activation products in a tape9 dictionary.
"""

ACTINIDE_FIELDS = ('sigma_alpha', 'sigma_f')
"""The cross section data keys for actinides & daughters in a tape9 dictionary.
"""

FISSION_PRODUCT_FIELDS = ('sigma_3n', 'sigma_p', 'TH232_fiss_yield',
                          'U233_fiss_yield',
                          'U235_fiss_yield', 'U238_fiss_yield',
                          'PU239_fiss_yield', 'PU241_fiss_yield',
                          'CM245_fiss_yield', 'CF249_fiss_yield')
"""The cross section & yield data keys for fission products in a tape9
   dictionary."""

# Table 4.2 in ORIGEN 2.2 manual
ORIGEN_TIME_UNITS = [None,              # No zero unit
                     1.0,               # seconds
                     60.0,              # minutes
                     3600.0,            # hours
                     86400.0,           # days
                     31556926.0,        # years...which are fuzzily defined.
                     np.inf,            # stable
                     31556926.0 * 1E3,  # ky
                     31556926.0 * 1E6,  # My
                     31556926.0 * 1E9,  # Gy
                     ]


def sec_to_time_unit(s):
    """Converts seconds to ORIGEN time and units.

    Parameters
    ----------
    s : float
        time in seconds

    Returns
    -------
    t : float
        time in units
    unit : int
        time unit that t is in. Represents index
        into ORIGEN_TIME_UNITS, which matches
        Table 4.2 in ORIGEN 2.2 manual.
    """
    for i, val in enumerate(ORIGEN_TIME_UNITS):
        if val is None:
            continue

        t = s / val
        unit = i

        if t != 0.0 and val == np.inf:
            # Origen spec for stable nuclides
            t = 0.0
            break
        elif 0.0 < t < 1.0:
            if i == 1:
                pass
            elif i == 7:
                unit -= 2
            else:
                unit -= 1
            t = s / ORIGEN_TIME_UNITS[unit]
            break

    return t, unit


# Regex helpers
_data_format = "\d+\.\d*[EeDd]?[ +-]?\d+"


###################################
### ORIGEN Input Deck Functions ###
###################################

def write_tape4(mat, outfile="TAPE4.INP"):
    """Writes a TAPE4.INP ORIGEN input file for a material.

    Parameters
    ----------
    mat : Material
        Material with mass weights in units of grams.
    outfile : str or file handler, optional
        Path to tape4 file or file-like object.
    """
    lower_z = mat[:'AC']
    upper_z = mat['AC':]

    lower_lines = ["1 {0} {1:.10E}   0 0   0 0   0 0".format(nucname.zzaaam(nuc), mass) \
                   for nuc, mass in lower_z.mult_by_mass().items()]
    upper_lines = ["2 {0} {1:.10E}   0 0   0 0   0 0".format(nucname.zzaaam(nuc), mass) \
                   for nuc, mass in upper_z.mult_by_mass().items()]
    lines = lower_lines + upper_lines + ["0 0 0 0\n"]

    tape4 = "\n".join(lines)

    # Write to the file
    opened_here = False
    if isinstance(outfile, basestring):
        outfile = open(outfile, 'w')
        opened_here = True

    outfile.write(tape4)

    if opened_here:
        outfile.close()


_tape5_irradiation_template = """\
  -1
  -1
  -1
  CUT     5 {CUT_OFF} -1
  RDA     Make sure thet the library identifier numbers match those in the TAPE9.INP file
  LIB     0 {DECAY_NLB1} {DECAY_NLB2} {DECAY_NLB3} {XSFPY_NLB1} {XSFPY_NLB2} {XSFPY_NLB3} 9 3 0 4 0
  OPTL    {optl}
  OPTA    {opta}
  OPTF    {optf}
  INP     1 -1  0  -1  4  4
  RDA     All irradiation (IRF and IRP) cards must be between burnup (BUP) cards.
  BUP
  {irr_type}     {irr_time}  {irr_value}   1   2   4  2
  BUP
  OUT     2  1 1 0
  END
"""

_tape5_decay_template = """\
  -1
  -1
  -1
  CUT     5 {CUT_OFF} -1
  RDA     Make sure thet the library identifier numbers match those in the TAPE9.INP file
  LIB     0 {DECAY_NLB1} {DECAY_NLB2} {DECAY_NLB3} {XSFPY_NLB1} {XSFPY_NLB2} {XSFPY_NLB3} 9 3 0 4 0
  OPTL    {optl}
  OPTA    {opta}
  OPTF    {optf}
  INP     1 -1  0  -1  4  4
  RDA     All irradiation (IRF and IRP) cards must be between burnup (BUP) cards.
  BUP
  DEC     {dec_time}  1   2   4  2
  BUP
  OUT     2  1 1 0
  END
"""

_nes_table = np.zeros((2, 2, 2), dtype=int)
_nes_table[True, True, True]    = 1
_nes_table[True, True, False]   = 2
_nes_table[True, False, True]   = 3
_nes_table[False, True, True]   = 4
_nes_table[True, False, False]  = 5
_nes_table[False, True, False]  = 6
_nes_table[False, False, True]  = 7
_nes_table[False, False, False] = 8


def _out_table_string(out_table_nes, out_table_num):
    """Makes a string output table line from relevant information."""
    if out_table_num is None:
        arr = np.ones(24, dtype=int)
        s = np.array2string(arr)[1:-1]
        return s

    arr = 8 * np.ones(24, dtype=int)
    idx = np.array(out_table_num, dtype=int) - 1

    arr[idx] = _nes_table[tuple(out_table_nes)]

    s = np.array2string(arr)[1:-1]

    return s


def write_tape5_irradiation(irr_type, irr_time, irr_value,
                            outfile="TAPE5.INP",
                            decay_nlb=(1, 2, 3),
                            xsfpy_nlb=(204, 205, 206),
                            cut_off=1E-10,
                            out_table_nes=(False, False, True),
                            out_table_laf=(True,  True,  True),
                            out_table_num=None):
    """Writes an irradiation TAPE5 file.

    Parameters
    ----------
    irr_type : str
        Flag that determines whether this is a constant power "IRP"
        irradiation or a constant flux "IRF" irradiation calculation.
    irr_time : float
        Irradiation time duration in days.
    irr_value : float
        Magnitude of the irradiation. If irr_type = "IRP", then
        this is a power.  If irr_type = "IRF", then this is a flux.
    outfile : str or file-like object
        Path or file to write the tape5 to.
    decay_nlb : length 3 sequence
        Three tuple of library numbers from the tape9 file decay data, eg (1, 2, 3).
    xsfpy_nlb : length 3 sequence
        Three tuple of library numbers from the tape9 file for cross section and fission
        product yields, eg (204, 205, 206).
    cut_off : float, optional
        Cut-off concentration, below which results are not recorded.
    out_table_nes :  length 3 sequence of bools, optional
        Specifies which type of output tables should be printed by ORIGEN.  The fields
        represent (Nuclide, Element, Summary).  The default value of (False, False, True)
        only prints the summary tables.
    out_table_laf :  length 3 sequence of bools, optional
        Specifies whether to print the activation products (l), actinides (a), and
        fission products (f).  By default all three are printed.
    out_table_num : sequence of ints or None
        Specifies which tables, by number, to print according to the rules given by
        out_table_nes and out_table_laf.  For example the list [10, 5] would print
        tables 5 and 10.  There are 24 tables available. If None, then all tables
        are printed.

    Warnings
    --------
    If ``irr_value`` is ``NaN`` or ``inf``, ORIGEN will still run without
    complaint, but the TAPE6.OUT file will only contain headers and no data.
    """
    if irr_type not in ["IRP", "IRF"]:
        raise TypeError("Irradiation type must be either 'IRP' or 'IRF'.")

    if np.isnan(irr_value):
        raise ValueError("Irradiation value is NaN.")

    if np.isinf(irr_value):
        raise ValueError("Irradiation value is infinite.")

    # Make template fill-value dictionary
    tape5_kw = {
        'CUT_OFF': "{0:.3E}".format(cut_off),
        'DECAY_NLB1': decay_nlb[0],
        'DECAY_NLB2': decay_nlb[1],
        'DECAY_NLB3': decay_nlb[2],
        'XSFPY_NLB1': xsfpy_nlb[0],
        'XSFPY_NLB2': xsfpy_nlb[1],
        'XSFPY_NLB3': xsfpy_nlb[2],
        'irr_type': irr_type,
        'irr_time': '{0:.10E}'.format(irr_time),
        'irr_value': '{0:.10E}'.format(irr_value),
        }

    no_print_string = np.array2string(8 * np.ones(24, dtype=int))[1:-1]

    # Activation Product Print String
    if out_table_laf[0]:
        tape5_kw['optl'] = _out_table_string(out_table_nes, out_table_num)
    else:
        tape5_kw['optl'] = no_print_string

    # Actinide Print String
    if out_table_laf[1]:
        tape5_kw['opta'] = _out_table_string(out_table_nes, out_table_num)
    else:
        tape5_kw['opta'] = no_print_string

    # Fission Product Print String
    if out_table_laf[2]:
        tape5_kw['optf'] = _out_table_string(out_table_nes, out_table_num)
    else:
        tape5_kw['optf'] = no_print_string

    # Fill the template and write it to a file
    tape5 = _tape5_irradiation_template.format(**tape5_kw)

    opened_here = False
    if isinstance(outfile, basestring):
        outfile = open(outfile, 'w')
        opened_here = True

    outfile.write(tape5)

    if opened_here:
        outfile.close()


def write_tape5_decay(dec_time,
                      outfile="TAPE5.INP",
                      decay_nlb=(1, 2, 3),
                      xsfpy_nlb=(204, 205, 206),
                      cut_off=1E-10,
                      out_table_nes=(False, False, True),
                      out_table_laf=(True,  True,  True),
                      out_table_num=None):
    """Writes a decay TAPE5 file.

    Parameters
    ----------
    dec_time : float
        Decay time duration in days.
    outfile : str or file-like object
        Path or file to write the tape5 to.
    decay_nlb : length 3 sequence
        Three tuple of library numbers from the tape9 file decay data, eg (1, 2, 3).
    xsfpy_nlb : length 3 sequence
        Three tuple of library numbers from the tape9 file for cross section and fission
        product yields, eg (204, 205, 206).
    cut_off : float, optional
        Cut-off concentration, below which reults are not recorded.
    out_table_nes :  length 3 sequence of bools, optional
        Specifies which type of output tables should be printed by ORIGEN.  The fields
        represent (Nuclide, Element, Summary).  The default value of (False, False, True)
        only prints the summary tables.
    out_table_laf :  length 3 sequence of bools, optional
        Specifies whether to print the activation products (l), actinides (a), and
        fission products (f).  By default all three are printed.
    out_table_num : sequence of ints or None
        Specifies which tables, by number, to print according to the rules given by
        out_table_nes and out_table_laf.  For example the list [10, 5] would print
        tables 5 and 10.  There are 24 tables available. If None, then all tables
        are printed.
    """
    # Make template fill-value dictionary
    tape5_kw = {
        'CUT_OFF': "{0:.3E}".format(cut_off),
        'DECAY_NLB1': decay_nlb[0],
        'DECAY_NLB2': decay_nlb[1],
        'DECAY_NLB3': decay_nlb[2],
        'XSFPY_NLB1': xsfpy_nlb[0],
        'XSFPY_NLB2': xsfpy_nlb[1],
        'XSFPY_NLB3': xsfpy_nlb[2],
        'dec_time': '{0:.10E}'.format(dec_time),
        }

    no_print_string = np.array2string(8 * np.ones(24, dtype=int))[1:-1]

    # Activation Product Print String
    if out_table_laf[0]:
        tape5_kw['optl'] = _out_table_string(out_table_nes, out_table_num)
    else:
        tape5_kw['optl'] = no_print_string

    # Actinide Print String
    if out_table_laf[1]:
        tape5_kw['opta'] = _out_table_string(out_table_nes, out_table_num)
    else:
        tape5_kw['opta'] = no_print_string

    # Fission Product Print String
    if out_table_laf[2]:
        tape5_kw['optf'] = _out_table_string(out_table_nes, out_table_num)
    else:
        tape5_kw['optf'] = no_print_string


    # Fill the template and write it to a file
    tape5 = _tape5_decay_template.format(**tape5_kw)

    opened_here = False
    if isinstance(outfile, basestring):
        outfile = open(outfile, 'w')
        opened_here = True

    outfile.write(tape5)

    if opened_here:
        outfile.close()

#
# Tape6 functions
#

_rx_bu_data_line = re.compile(' (TIME, SEC|NEUT. FLUX|SP POW,MW|BURNUP,MWD|K INFINITY|NEUT PRODN|NEUT DESTN|TOT BURNUP|AVG N FLUX|AVG SP POW) (.*)')

_rx_bu_key_map = {
    "TIME, SEC":  "time_sec",
    "NEUT. FLUX": "flux",
    "SP POW,MW":  "specific_power_MW",
    "BURNUP,MWD": "burnup_MWD",
    "K INFINITY": "k_inf",
    "NEUT PRODN": "neutron_production_rate",
    "NEUT DESTN": "neutron_destruction_rate",
    "TOT BURNUP": "total_burnup",
    "AVG N FLUX": "average_flux",
    "AVG SP POW": "average_specific_power",
    }

_table_header_line = re.compile("[ 0]\s+(\d+) (NUCLIDE|ELEMENT|SUMMARY) TABLE:([ A-Z]+),([ 0-9A-Za-z*]+)")

_table_header_alpha_line = re.compile("[ 0]\s+(\d+) (NUCLIDE|ELEMENT|SUMMARY) TABLE:\s*(ALPHA RADIOACTIVITY)\s+([ 0-9A-Za-z*]+)")

_nuclide_line = re.compile(" ([ A-Z][A-Z][ \d][ \d]\d[ M])\s+(.*)")

_element_line = re.compile(" ([ A-Z][A-Z])   \s+(.*)")

_species_group_line = re.compile('[ +]\s+(ACTIVATION PRODUCTS|ACTINIDES[+]DAUGHTERS|FISSION PRODUCTS)')

_group_key_map = {
    'ACTIVATION PRODUCTS': 'activation_products',
    'ACTINIDES+DAUGHTERS': 'actinides',
    'FISSION PRODUCTS': 'fission_products',
    }

_alpha_n_header_line = re.compile('\s*(\(ALPHA,N\) NEUTRON SOURCE), (NEUTRONS/SEC)')

_spont_fiss_header_line = re.compile('\s*(SPONTANEOUS FISSION NEUTRON SOURCE), (NEUTRONS/SEC)')

_n_source_key_map = {
    '(ALPHA,N) NEUTRON SOURCE': 'alpha_neutron_source',
    'SPONTANEOUS FISSION NEUTRON SOURCE': 'spont_fiss_neutron_source',
    }

_photon_spec_header_line = re.compile("\s+PHOTON SPECTRUM FOR(.*)")


def parse_tape6(tape6="TAPE6.OUT"):
    """Parses an ORIGEN 2.2 TAPE6.OUT file.

    Parameters
    ----------
    tape6 : str or file-like object
        Path or file to read the tape6 file from.

    Returns
    -------
    results : dict
        Dictionary of parsed values.

    Warnings
    --------
    This method currently only functions to extract neutronic data from TAPE6
    files.  It does not yet parse out photonic data.  If you would like to see
    this feature added, please contact the developers.

    Notes
    -----
    The results dictionary that is returned is highly structured and generally
    matches the layout of the TAPE6 file.  Data is stored as 1d numpy float arrays
    which (if the TAPE6 is well-formed) will all be of the same length and match
    the time vector.  The possible layout of results is as follows::

      |- 'time_sec': time per index in [seconds]
      |- 'flux': neutron flux at this time [n/cm^2/s]
      |- 'specific_power_MW': recator specific power at this time [MW]
      |- 'burnup_MWD': reactor burnup since last time step [MWd/input mass [g] from TAPE4]
      |- 'k_inf': infinite multiplication factor [unitless]
      |- 'neutron_production_rate': Total reactor neutron production rate [n/s]
      |- 'neutron_destruction_rate: Total reactor neutron destruction rate [n/s]
      |- 'total_burnup': Cummulative burnup over all time [MWd/input mass [g] from TAPE4]
      |- 'average_flux': average neutron flux over preceeding time interval [n/cm^2/s]
      |- 'average_specific_power: recator specific power over preceeding time interval [MW]
      |- 'materials': list of Materials of same length as 'time_sec', only present if
      |               'table_3' or 'table_5' exist and have 'nuclide' output.
      |- 'alpha_neutron_source': dict
      |                          |- 'title': str
      |                          |- 'units': str
      |                          |- nuclide or element str: (alpha, n) neutron source [n/s]
      |- 'spont_fiss_neutron_source': dict
      |                          |- 'title': str
      |                          |- 'units': str
      |                          |- nuclide or element str: spontaneous fission neutron source [n/s]
      |- 'table_{n}': dict
      |               |- 'nuclide': dict
      |               |             |- 'title': str
      |               |             |- 'units': str
      |               |             |- 'activation_products': dict of (nuc-zzaaam, data) pairs
      |               |             |- 'actinides': dict of (nuc-zzaaam, data) pairs
      |               |             |- 'fission_products': dict of (nuc-zzaaam, data) pairs
      |               |- 'element': dict
      |               |             |- 'title': str
      |               |             |- 'units': str
      |               |             |- 'activation_products': dict of (elem str, data) pairs
      |               |             |- 'actinides': dict of (elem str, data) pairs
      |               |             |- 'fission_products': dict of (elem str, data) pairs
      |               |- 'summary': dict
      |               |             |- 'title': str
      |               |             |- 'units': str
      |               |             |- 'activation_products': dict of (elem or nuc str, data) pairs
      |               |             |- 'actinides': dict of (elem or nuc str, data) pairs
      |               |             |- 'fission_products': dict of (elem or nuc str, data) pairs

    """
    # Read the TAPE6 file
    opened_here = False
    if isinstance(tape6, basestring):
        tape6 = open(tape6, 'r')
        opened_here = True

    lines = tape6.readlines()

    if opened_here:
        tape6.close()

    # Prep to parse the file
    results = {}

    # Defaults
    table_key = None
    table_type = None
    table_group = None

    # Read in the file line-by-line
    for i, line in enumerate(lines):
        # Get reactivity and burnup data
        m = _rx_bu_data_line.match(line)
        if m is not None:
            key, data = m.groups()
            new_key = _rx_bu_key_map[key]
            arr_data = np.array(data.split(), dtype=float)
            curr_data = results.get(new_key, [])
            results[new_key] = np.append(curr_data, arr_data)
            continue

        # Get table spcies group
        m = _species_group_line.match(line)
        if m is not None:
            table_group = _group_key_map[m.group(1)]
            continue

        # Get table header info
        m = _table_header_line.match(line) or _table_header_alpha_line.match(line)
        if m is not None:
            tnum, ttype, ttitle, tunits = m.groups()

            table_key = "table_{0}".format(tnum)
            if table_key not in results:
                results[table_key] = {}

            table_type = ttype.lower()
            if table_type not in results[table_key]:
                results[table_key][table_type] = {}

            results[table_key][table_type]["title"] = ttitle.strip().lower()
            results[table_key][table_type]["units"] = tunits.strip().lower()
            if table_group not in results[table_key][table_type]:
                results[table_key][table_type][table_group] = {}
            continue

        # Grab nuclide data lines
        m = _nuclide_line.match(line)
        if (m is not None) and (table_key is not None):
            nuc, data = m.groups()
            nuc_name = nuc.replace(' ', '')

            # Don't know WTF element 'SF' is suppossed to be!
            # (Spent fuel, spontaneous fission)
            if nuc_name == 'SF250':
                continue

            nuc_zz = nucname.zzaaam(nuc_name)
            nuc_key = nuc_zz if table_type == 'nuclide' else nuc_name
            nuc_data = np.array(data.split(), dtype=float)

            if table_key.startswith('table_'):
                curr_data = results[table_key][table_type][table_group].get(nuc_key, [])
                results[table_key][table_type][table_group][nuc_key] = np.append(curr_data, nuc_data)
            else:
                curr_data = results[table_key].get(nuc_key, [])
                results[table_key][nuc_key] = np.append(curr_data, nuc_data)
            continue

        # Grab element data line
        m = _element_line.match(line)
        if (m is not None) and (table_key is not None):
            elem, data = m.groups()
            elem = elem.replace(' ', '')

            # Still don't know WTF element 'SF' is suppossed to be!
            # (Spent fuel, spontaneous fission)
            if elem == 'SF':
                continue

            elem_data = np.array(data.split(), dtype=float)

            if table_key.startswith('table_'):
                curr_data = results[table_key][table_type][table_group].get(elem, [])
                results[table_key][table_type][table_group][elem] = np.append(curr_data, elem_data)
            else:
                curr_data = results[table_key].get(elem, [])
                results[table_key][elem] = np.append(curr_data, elem_data)
            continue

        # Grab (alpha, n) and spontaneous fission headers
        m = _alpha_n_header_line.match(line) or _spont_fiss_header_line.match(line)
        if m is not None:
            ttitle, tunits = m.groups()

            table_key = _n_source_key_map[ttitle]
            if table_key not in results:
                results[table_key] = {}

            table_type = None
            table_group = None

            results[table_key]["title"] = ttitle.strip().lower()
            results[table_key]["units"] = tunits.strip().lower()
            continue

        # Photon spectra parsing is not yet supported
        m = _photon_spec_header_line.match(line)
        if m is not None:
            table_key = None
            table_type = None
            table_group = None

    # Done with parsing, try to convert to material
    tbl = None
    if ('table_5' in results) and ('nuclide' in results['table_5']):
        tbl = 'table_5'
        mat_gen = Material
    elif ('table_3' in results) and ('nuclide' in results['table_3']):
        tbl = 'table_3'
        mat_gen = from_atom_frac

    if tbl is not None:
        T = len(results['time_sec'])
        mats = [Material() for t in range(T)]

        for grp in _group_key_map.values():
            if grp in results[tbl]['nuclide']:
                mats = [m + mat_gen(dict([(nuc, arr[i]) for nuc, arr in \
                                            results[tbl]['nuclide'][grp].items()])) \
                                            for i, m in enumerate(mats)]

        results['materials'] = mats

    return results


#
# Tape9 functions
#

title_card_re = re.compile("(\d+)\s+(\S.*)")

# Decay library regex
decay_card1_re = re.compile("(\d+)\s+(\d{{5,7}})\s+(\d)\s+({num})\s+({num})\s+({num})\s+({num})\s+({num})\s+({num})".format(num=_data_format))
decay_card2_re = re.compile("(\d+)\s+({num})\s+({num})\s+({num})\s+({num})\s+({num})\s+({num})".format(num=_data_format))

# Cross section and fission product yeild library regex
xsfpy_card1_re = re.compile("(\d+)\s+(\d{{5,7}})\s+({num})\s+({num})\s+({num})\s+({num})\s+({num})\s+({num})\s+([+-]?{num})".format(num=_data_format))
xsfpy_card2_re = re.compile("(\d+)\s+({num})\s+({num})\s+({num})\s+({num})\s+({num})\s+({num})\s+({num})\s+({num})".format(num=_data_format))


def _parse_tape9_decay(deck):
    pdeck = {'_type': 'decay'}
    pdeck['title'] = title_card_re.match(deck[0]).group(2).strip()

    # Parse the cards into a structured arrau
    cards = [m.groups()[1:] + n.groups()[1:] for m, n in
             zip(map(decay_card1_re.match, deck[1::2]),
                 map(decay_card2_re.match, deck[2::2]))]
    cards = [tuple(d.replace(' ', '') for d in card) for card in cards]
    cards = np.array(cards, dtype='i4,i4' + ',f8'*12)
    pdeck['_cards'] = cards

    # Add the first cards
    pdeck['half_life'] = dict([(nuc, ORIGEN_TIME_UNITS[unit]*(val or 1.0)) for nuc, unit, val in cards[['f0', 'f1', 'f2']] ])
    pdeck['frac_beta_minus_x'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f3']] ])
    pdeck['frac_beta_plus_or_electron_capture'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f4']] ])
    pdeck['frac_beta_plus_or_electron_capture_x'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f5']] ])
    pdeck['frac_alpha'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f6']] ])
    pdeck['frac_isomeric_transition'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f7']] ])

    # Add the second cards
    pdeck['frac_spont_fiss'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f8']] ])
    pdeck['frac_beta_n'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f9']] ])
    pdeck['recoverable_energy'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f10']] ])
    pdeck['frac_natural_abund'] = dict([(nuc, val*0.01) for nuc, val in cards[['f0', 'f11']] ])
    pdeck['inhilation_concentration'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f12']] ])
    pdeck['ingestion_concentration'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f13']] ])

    return pdeck


def _parse_tape9_xsfpy(deck):
    pdeck = {'_type': 'xsfpy'}
    pdeck['title'] = title_card_re.match(deck[0]).group(2).strip()

    # Pasre the deck
    cards = []
    no_fpy = ('0.0', ) * 8

    i = 1
    deck_size = len(deck)
    while (i < deck_size):
        first_card = xsfpy_card1_re.match(deck[i]).groups()[1:]
        if 0.0 < float(first_card[-1]):
            i += 1
            second_card = xsfpy_card2_re.match(deck[i]).groups()[1:]
        else:
            second_card = no_fpy
        cards.append(first_card + second_card)
        i += 1

    cards = [tuple(d.replace(' ', '') for d in card) for card in cards]
    cards = np.array(cards, dtype='i4' + ',f8'*15)
    pdeck['_cards'] = cards

    # Try to determine subtype
    if (0.0 < cards['f7']).any():
        subtype = 'fission_products'
    elif 890000 < cards['f0'].mean():
        subtype = 'actinides'
    else:
        subtype = 'activation_products'
    pdeck['_subtype'] = subtype

    # Parse first cards
    pdeck['sigma_gamma'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f1']] ])
    pdeck['sigma_2n'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f2']] ])

    f3_keys = {'fission_products': 'sigma_alpha', 'actinides': 'sigma_3n', 'activation_products': 'sigma_alpha'}
    pdeck[f3_keys[subtype]] = dict([(nuc, val) for nuc, val in cards[['f0', 'f3']] ])

    f4_keys = {'fission_products': 'sigma_p', 'actinides': 'sigma_f', 'activation_products': 'sigma_p'}
    pdeck[f4_keys[subtype]] = dict([(nuc, val) for nuc, val in cards[['f0', 'f4']] ])

    pdeck['sigma_gamma_x'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f5']] ])
    pdeck['sigma_2n_x'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f6']] ])

    pdeck['fiss_yields_present'] = dict([(nuc, 0.0 < val) for nuc, val in cards[['f0', 'f7']] ])

    # parse second cards if of correct subtype
    if subtype == 'fission_products':
        pdeck['TH232_fiss_yield'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f8']] ])
        pdeck['U233_fiss_yield'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f9']] ])
        pdeck['U235_fiss_yield'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f10']] ])
        pdeck['U238_fiss_yield'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f11']] ])
        pdeck['PU239_fiss_yield'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f12']] ])
        pdeck['PU241_fiss_yield'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f13']] ])
        pdeck['CM245_fiss_yield'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f14']] ])
        pdeck['CF249_fiss_yield'] = dict([(nuc, val) for nuc, val in cards[['f0', 'f15']] ])

    return pdeck


def parse_tape9(tape9="TAPE9.INP"):
    """Parses an ORIGEN 2.2 TAPE9 file and returns the data as a dictionary of
    nuclide dictionaries.

    Parameters
    ----------
    tape9 : str or file-like object, optional
        Path to the tape9 file.

    Returns
    -------
    parsed : dict
        A dictionary of the data from the TAPE9 file.

    Notes
    -----
    The TAPE9 format is highly structured. Therefore the in-memory representation contains a
    non-trivial amount of nesting.  At the top level, the dictionary keys are the library
    numbers::

        tape9
          |- keys : deck number (1, 2, 3, 241, ...)
          |- values : sub-dictionaries for each deck.

    Each deck contains keys which vary by deck type (and subtype).  All dictionary-typed data
    maps zzaaam-nuclide integers to the appropriate value::

        all decks
          |- '_type' : str in ['decay', 'xsfpy'] # photon libs not yet supported
          |- '_subtype' : str for 'xsfpy' in ['activation_products', 'actinides', 'fission_products']
          |- 'title' : str, deck name
          |- '_cards' : optional, numpy structrued array of deck data

        decay decks
          |- 'half_life' : float-valued dict [seconds]
          |- 'frac_beta_minus_x' : float-valued dict [fraction of decays via beta minus
          |                        which leave an excited nucleus]
          |- 'frac_beta_plus_or_electron_capture' : float-valued dict [fraction of decays
          |                                         via positron emission or electron capture]
          |- 'frac_beta_plus_or_electron_capture_x' : float-valued dict [fraction of decays
          |                                           via positron emission or electron capture
          |                                           which leave an excited nucleus]
          |- 'frac_alpha' : float-valued dict [fraction of decays via alpha emission]
          |- 'frac_isomeric_transition' : float-valued dict [fraction of decays from an excitied
          |                               state to the ground state]
          |- 'frac_spont_fiss' : float-valued dict [fraction of decays via spontateous fission]
          |- 'frac_beta_n' : float-valued dict [fraction of decays via beta plus a neutron]
          |- 'recoverable_energy' : float-valued dict, Total recoverable energy [MeV / decay]
          |- 'frac_natural_abund' : float-valued dict, natrual occuring abundance [atom fraction]
          |- 'inhilation_concentration' : float-valued dict, continuous inhilation [RCG]
          |- 'ingestion_concentration' : float-valued dict, continuous ingestion [RCG]

        cross section and fission product yield decks
          |- 'sigma_gamma' : float-valued dict, (n, gamma) cross section [barns]
          |- 'sigma_2n' : float-valued dict, (n, 2n) cross section [barns]
          |- 'sigma_gamma_x' : float-valued dict, (n, gamma *) cross section [barns]
          |- 'sigma_2n_x' : float-valued dict, (n, 2n *) cross section [barns]
          |- 'fiss_yields_present' : bool-valued dict, Whether fission product yields are
                                     included for this nuclide.

        activation product cross section decks
          |- '_subtype' : 'activation_products'
          |- 'sigma_3n' : float-valued dict, (n, 3n) cross section [barns]
          |- 'sigma_p' : float-valued dict, (n, proton) cross section [barns]

        actinide cross section decks
          |- '_subtype' : 'actinides'
          |- 'sigma_alpha' : float-valued dict, (n, alpha) cross section [barns]
          |- 'sigma_f' : float-valued dict, (n, fission) cross section [barns]

        fission product cross section and yield decks
          |- '_subtype' : 'fission_products'
          |- 'sigma_3n' : float-valued dict, (n, 3n) cross section [barns]
          |- 'sigma_p' : float-valued dict, (n, proton) cross section [barns]
          |- 'TH232_fiss_yield' : float-valued dict, yield from Th-232 fission [frac]
          |- 'U233_fiss_yield' : float-valued dict, yield from U-233 fission [frac]
          |- 'U235_fiss_yield' : float-valued dict, yield from U-235 fission [frac]
          |- 'U238_fiss_yield' : float-valued dict, yield from U-238 fission [frac]
          |- 'PU239_fiss_yield' : float-valued dict, yield from Pu-239 fission [frac]
          |- 'PU241_fiss_yield' : float-valued dict, yield from Pu-241 fission [frac]
          |- 'CM245_fiss_yield' : float-valued dict, yield from Cm-245 fission [frac]
          |- 'CF249_fiss_yield' : float-valued dict, yield from Cf-249 fission [frac]
    """
    # Read and strip lines
    opened_here = False
    if isinstance(tape9, basestring):
        tape9 = open(tape9, 'r')
        opened_here = True

    tape9_lines = [line.strip() for line in tape9]

    if opened_here:
        tape9.close()

    # Split lines into various decks.
    decks = []
    while 0 < tape9_lines.count('-1'):
        n = tape9_lines.index('-1')
        decks.append(tape9_lines[:n])
        tape9_lines = tape9_lines[n+1:]

    # parse the individual decks.
    parsed = {}
    for deck in decks:
        # Decay deck
        m = decay_card1_re.match(deck[1])
        if m is not None:
            deck_num = int(m.group(1))
            parsed[deck_num] = _parse_tape9_decay(deck)
            continue

        # Cross section deck
        m = xsfpy_card1_re.match(deck[1])
        if m is not None:
            deck_num = int(m.group(1))
            parsed[deck_num] = _parse_tape9_xsfpy(deck)
            continue

    return parsed


def loads_tape9(tape9):
    """Parses a string that represents an ORIGEN 2.2 TAPE9 file and returns the data
    as a dictionary of nuclide dictionaries. See ``parse_tape9()`` for more details.

    Parameters
    ----------
    tape9 : str
        String represetation of the TAPE9 file.

    Returns
    -------
    parsed : dict
        A dictionary of the data from the TAPE9 file.
    """
    if isinstance(tape9, unicode):
        t9 = StringIO(tape9)
    else:
        t9 = StringIO(tape9.decode())
    parsed = parse_tape9(t9)
    return parsed


def merge_tape9(tape9s):
    """Merges a sequence of full or partial TAPE9s into a single tape9 dictionary.
    Data from the first tape9 has precednce over the second, the second over the
    third, etc.

    Parameters
    ----------
    tape9s : list or tuple
        A sequence of full or partial tape9 dictionaries.  See parse_tape9()
        for more information on the structure of these dictionaires.

    Returns
    -------
    tape9 : dictionary
        A tape9 file which is the merger of the all of the tape9s.  See
        parse_tape9() for more information on the structure of this dictionary.
    """
    tape9 = {}

    for t9 in tape9s[::-1]:
        for nlb, deck in t9.items():
            if nlb in tape9:
                # Make sure the decks are of the same type
                assert tape9[nlb]['_type'] == deck['_type']
                if ('_subtype' in tape9[nlb]) and ('_subtype' in deck):
                    assert tape9[nlb]['_subtype'] == deck['_subtype']

                # _cards has been invalidated... remove it
                tape9[nlb].pop('_cards', None)

                # Update all of the keys, except _cards
                for key, value in deck.items():
                    if key in tape9[nlb] and hasattr(value, 'keys'):
                        tape9[nlb][key].update(value)
                    elif key != '_cards':
                        tape9[nlb][key] = deepcopy(value)
            else:
                # New library number, make a copy
                tape9[nlb] = deepcopy(deck)

    return tape9


def _double_get(dict, key1, key2, default=0.0):
    if key1 in dict:
        return dict[key1].get(key2, default)
    else:
        return default


_deck_title_fmt = "{nlb:>4}    {title:^72}\n"

_decay_card_fmt = ("{nlb:>4}{nuc:>8}  {unit}     {time:<9.{p}E} {fbx:<9.{p}E} "
                   "{fpec:<9.{p}E} {fpecx:<9.{p}E} {fa:<9.{p}E} {fit:<9.{p}E}\n"
                   "{nlb:>4}                {fsf:<9.{p}E} {fn:<9.{p}E} "
                   "{qrec:<9.{p}E} {abund:<9.{p}E} {arcg:<9.{p}E} "
                   "{wrcg:<9.{p}E}\n")

_xs_card_fmt = ("{nlb:>4}{nuc:>8} {sg:<9.{p}E} {s2n:<9.{p}E} {s3n_or_a:<9.{p}E} "
                "{sf_or_p:<9.{p}E} {sg_x:<9.{p}E} {s2n_x:<9.{p}E} "
                "{fpy_flag:>6.1F} \n")

_fpy_card_fmt = ("{nlb:>4}     {y1:<8.{p}E} {y2:<8.{p}E} {y3:<8.{p}E} "
                 "{y4:<8.{p}E} {y5:<8.{p}E} {y6:<8.{p}E} {y7:<8.{p}E} "
                 "{y8:<8.{p}E}\n")


def _decay_deck_2_str(nlb, deck, precision):
    # Get unique isotopes
    nucset = set([nuc for nuc in chain(*[v.keys() for k, v in deck.items() \
                  if hasattr(v, 'keys')]) ])
    nucset = sorted(map(int, nucset))

    s = ""
    for nuc in nucset:
        nuc_id = nucname.zzaaam_to_id(nuc)
        t, unit = sec_to_time_unit(_double_get(deck, 'half_life', nuc,
                                               data.half_life(nuc_id)))
        s += _decay_card_fmt.format(nlb=nlb,
                nuc=nuc,
                unit=unit,
                time=t,
                fbx=_double_get(deck, 'frac_beta_minus_x', nuc),
                fpec=_double_get(deck, 'frac_beta_plus_or_electron_capture', nuc),
                fpecx=_double_get(deck, 'frac_beta_plus_or_electron_capture_x', nuc),
                fa=_double_get(deck, 'frac_alpha', nuc),
                fit=_double_get(deck, 'frac_isomeric_transition', nuc),
                fsf=_double_get(deck, 'frac_spont_fiss', nuc),
                fn=_double_get(deck, 'frac_beta_n', nuc),
                qrec=_double_get(deck, 'recoverable_energy', nuc),
                abund=_double_get(deck, 'frac_natural_abund', nuc,
                                  data.natural_abund(nuc_id)),
                arcg=_double_get(deck, 'inhilation_concentration', nuc, 1.0),
                wrcg=_double_get(deck, 'ingestion_concentration', nuc, 1.0),
                p=precision,
                )
    return s


def _xs_deck_2_str(nlb, deck, precision):
    # Get unique isotopes
    nucset = set([nuc for nuc in chain(*[v.keys() for k, v in deck.items() if hasattr(v, 'keys')]) ])
    nucset = sorted(nucset)

    is_actinides = deck['_subtype'] == 'actinides'

    s = ""
    for nuc in nucset:
        fpy_flag = -1.0
        fpy_present = _double_get(deck, 'fiss_yields_present', nuc, False)
        if fpy_present:
            fpy_flag = 1.0

        s += _xs_card_fmt.format(nlb=nlb,
                                 nuc=nuc,
                                 sg=_double_get(deck, 'sigma_gamma', nuc),
                                 s2n=_double_get(deck, 'sigma_2n', nuc),
                                 s3n_or_a=_double_get(deck, 'sigma_alpha', nuc) if is_actinides else _double_get(deck, 'sigma_3n', nuc),
                                 sf_or_p=_double_get(deck, 'sigma_f', nuc) if is_actinides else _double_get(deck, 'sigma_p', nuc),
                                 sg_x=_double_get(deck, 'sigma_gamma_x', nuc),
                                 s2n_x=_double_get(deck, 'sigma_2n_x', nuc),
                                 fpy_flag=fpy_flag,
                                 p=precision,
                                 )
    return s


def _xsfpy_deck_2_str(nlb, deck, precision):
    # Get unique isotopes
    nucset = set([nuc for nuc in chain(*[v.keys() for k, v in deck.items() if hasattr(v, 'keys')]) ])
    nucset = sorted(nucset)
    s = ""
    fpy_precision = precision-1  # otherwise it doesn't fit in 80 chars
    for nuc in nucset:
        fpy_flag = 1.0
        s += _xs_card_fmt.format(nlb=nlb,
                                 nuc=nuc,
                                 sg=_double_get(deck, 'sigma_gamma', nuc),
                                 s2n=_double_get(deck, 'sigma_2n', nuc),
                                 s3n_or_a=_double_get(deck, 'sigma_3n', nuc),
                                 sf_or_p=_double_get(deck, 'sigma_p', nuc),
                                 sg_x=_double_get(deck, 'sigma_gamma_x', nuc),
                                 s2n_x=_double_get(deck, 'sigma_2n_x', nuc),
                                 fpy_flag=fpy_flag,
                                 p=precision,
                                 )
        s += _fpy_card_fmt.format(nlb=nlb,
                                  y1=_double_get(deck, 'TH232_fiss_yield', nuc),
                                  y2=_double_get(deck, 'U233_fiss_yield', nuc),
                                  y3=_double_get(deck, 'U235_fiss_yield', nuc),
                                  y4=_double_get(deck, 'U238_fiss_yield', nuc),
                                  y5=_double_get(deck, 'PU239_fiss_yield', nuc),
                                  y6=_double_get(deck, 'PU241_fiss_yield', nuc),
                                  y7=_double_get(deck, 'CM245_fiss_yield', nuc),
                                  y8=_double_get(deck, 'CF249_fiss_yield', nuc),
                                  p=fpy_precision,
                                  )
    return s


def _del_deck_nuc(deck, nuc):
    """removes a nucide from a deck completely."""
    for field, data in deck.items():
        if not isinstance(data, Mapping):
            continue
        if nuc not in data:
            continue
        del data[nuc]


def _filter_fpy(tape9):
    decay_nlb, xsfpy_nlb = nlbs(tape9)
    declib = tape9[decay_nlb[-1]]
    fpylib = tape9[xsfpy_nlb[-1]]
    nucset = set([nuc for nuc in chain(*[v.keys() for k, v in fpylib.items() \
                  if isinstance(v, Mapping)])])
    decnucs = set([nuc for nuc in chain(*[v.keys() for k, v in declib.items() \
                   if isinstance(v, Mapping)])])
    for nuc in nucset:
        fpy_present = _double_get(fpylib, 'fiss_yields_present', nuc, False)
        y1 = _double_get(fpylib, 'TH232_fiss_yield', nuc)
        y2 = _double_get(fpylib, 'U233_fiss_yield', nuc)
        y3 = _double_get(fpylib, 'U235_fiss_yield', nuc)
        y4 = _double_get(fpylib, 'U238_fiss_yield', nuc)
        y5 = _double_get(fpylib, 'PU239_fiss_yield', nuc)
        y6 = _double_get(fpylib, 'PU241_fiss_yield', nuc)
        y7 = _double_get(fpylib, 'CM245_fiss_yield', nuc)
        y8 = _double_get(fpylib, 'CF249_fiss_yield', nuc)
        fpy_present = fpy_present and any([y > 0.0 for y in [y1, y2, y3, y4,
                                                             y5, y6, y7, y8]])
        if fpy_present:
            continue
        _del_deck_nuc(fpylib, nuc)
        _del_deck_nuc(declib, nuc)


def _ensure_nucs_in_decay(tape9):
    decay_nlb, xsfpy_nlb = nlbs(tape9)
    for dn, xn in zip(decay_nlb, xsfpy_nlb):
        dlib = tape9[dn]
        dhl = tape9[dn]['half_life']
        xlib = tape9[xn]
        nucset = set([nuc for nuc in chain(*[v.keys() for k, v in xlib.items() \
                      if isinstance(v, Mapping)])])
        for nuc in nucset:
            if nuc not in dhl:
                dhl[nuc] = data.half_life(nucname.zzaaam(int(nuc)))


_DECK_2_STR_MAP = {
    ('decay', None): _decay_deck_2_str,
    ('xsfpy', 'activation_products'): _xs_deck_2_str,
    ('xsfpy', 'actinides'): _xs_deck_2_str,
    ('xsfpy', 'fission_products'): _xsfpy_deck_2_str,
    }


def write_tape9(tape9, outfile="TAPE9.INP", precision=3):
    """Writes an ORIGEN 2.2 TAPE9.INP file given a tape9 dictionary of values.

    Parameters
    ----------
    tape9 : dict
        A tape9 dictionary. See parse_tape9() for more information on the structure.
    outfile : str or file-like object, optional
        Path to the new tape9 file.
    precision :  int, optional
        The number of significant figures that all output data is given to beyond
        the decimal point.
    """
    t9 = ""
    _filter_fpy(tape9)
    _ensure_nucs_in_decay(tape9)
    for nlb, deck in tape9.items():
        t9 += _deck_title_fmt.format(nlb=nlb, title=deck['title'])
        t9 += _DECK_2_STR_MAP[deck['_type'], deck.get('_subtype', None)](nlb, deck, precision)
        t9 += "  -1\n"

    opened_here = False
    if isinstance(outfile, basestring):
        outfile = open(outfile, 'w')
        opened_here = True

    outfile.write(t9)

    if opened_here:
        outfile.close()

_fyp_present = {'activation_products': False,  'actinides': False,
                'fission_products': True, }

_xslib_computers = {
    'sigma_gamma': lambda nuc, xscache: xscache[nuc, 'gamma'][0],
    'sigma_2n': lambda nuc, xscache: xscache[nuc, 'z_2n'][0],
    'sigma_gamma_x': lambda nuc, xscache: xscache[nuc, 'gamma_1'][0] + \
                                          xscache[nuc, 'gamma_2'][0],
    'sigma_2n_x': lambda nuc, xscache: xscache[nuc, 'z_2n_1'][0] + \
                                       xscache[nuc, 'z_2n_2'][0],
    'sigma_3n': lambda nuc, xscache: xscache[nuc, 'z_3n'][0],
    'sigma_p': lambda nuc, xscache: xscache[nuc, 'p'][0],
    'sigma_alpha': lambda nuc, xscache: xscache[nuc, 'alpha'][0],
    'sigma_f': lambda nuc, xscache: xscache[nuc, 'fission'][0],
    'TH232_fiss_yield': lambda nuc, xscache: data.fpyield(902320000, nuc),
    'U233_fiss_yield': lambda nuc, xscache: data.fpyield(922330000, nuc),
    'U235_fiss_yield': lambda nuc, xscache: data.fpyield(922350000, nuc),
    'U238_fiss_yield': lambda nuc, xscache: data.fpyield(922380000, nuc),
    'PU239_fiss_yield': lambda nuc, xscache: data.fpyield(942390000, nuc),
    'PU241_fiss_yield': lambda nuc, xscache: data.fpyield(942410000, nuc),
    'CM245_fiss_yield': lambda nuc, xscache: data.fpyield(962450000, nuc),
    'CF249_fiss_yield': lambda nuc, xscache: data.fpyield(982490000, nuc),
    }


def _compute_xslib(nuc, key, lib, xscache):
    for field, data in lib.items():
        if field.startswith('_'):
            continue
        elif field == 'fiss_yields_present':
            data[key] = _fyp_present[lib['_subtype']]
            continue
        elif field == 'title':
            continue
        data[key] = _xslib_computers[field](nuc, xscache)


def xslibs(nucs=NUCS, xscache=None, nlb=(201, 202, 203), verbose=False):
    """Generates a TAPE9 dictionary of cross section & fission product yield data
    for a set of nuclides.

    Parameters
    ----------
    nucs : iterable of ints, optional
        Set of nuclides in id form
    xscache : XSCache, optional
        A cross section cache to get cross section data. If None, uses default.
    nlb : length-3 sequence of ints
        Library numbers for activation products, actinides & daugthers, and
        fission products respectively.
    verbose : bool, optional
        Flag to print status as we go.

    Returns
    -------
    t9 : dict
        The data needed for a TAPE9 file.
    """
    if xscache is None:
        xscache = cache.xs_cache
    old_flux = xscache.get('phi_g', None)
    old_group_struct = xscache.get('E_g', None)
    if old_group_struct is None:
        xscache['E_g'] = [10.0, 1e-7]
    elif len(old_group_struct) == 2:
        pass
    else:
        xscache['E_g'] = [old_group_struct[0], old_group_struct[-1]]
    nucs = sorted(nucs)

    # setup tape9
    t9 = {nlb[0]: {'_type': 'xsfpy', '_subtype': 'activation_products',
                   'title': 'PyNE Cross Section Data for Activation Products'},
          nlb[1]: {'_type': 'xsfpy', '_subtype': 'actinides',
                   'title': 'PyNE Cross Section Data for Actinides & Daughters'},
          nlb[2]: {'_type': 'xsfpy', '_subtype': 'fission_products',
                   'title': 'PyNE Cross Section Data for Fission Products'},
          }
    for n, lib in t9.items():
        for field in XSFPY_FIELDS:
            lib[field] = {}
    for field in ACTIVATION_PRODUCT_FIELDS:
        t9[nlb[0]][field] = {}
    for field in ACTINIDE_FIELDS:
        t9[nlb[1]][field] = {}
    for field in FISSION_PRODUCT_FIELDS:
        t9[nlb[2]][field] = {}

    # fill with data
    for nuc in nucs:
        if verbose:
            print('computing {0}'.format(nucname.name(nuc)))
        key = nucname.zzaaam(nuc)
        if nuc in ACTIVATION_PRODUCT_NUCS:
            _compute_xslib(nuc, key, t9[nlb[0]], xscache)
        if nuc in ACTINIDE_AND_DAUGHTER_NUCS:
            _compute_xslib(nuc, key, t9[nlb[1]], xscache)
        if nuc in FISSION_PRODUCT_NUCS:
            _compute_xslib(nuc, key, t9[nlb[2]], xscache)
    xscache['E_g'] = old_group_struct
    xscache['phi_g'] = old_flux
    return t9


def nlbs(t9):
    """Finds the library number tuples in a tape9 dictionary.

    Parameters
    ----------
    t9 : dict
        TAPE9 dictionary.

    Returns
    -------
    decay_nlb : 3-tuple
        Tuple of decay library numbers.
    xsfpy_nlb : 3-tuple
        Tuple of cross section & fission product library numbers.
    """
    decay_nlb = []
    xsfpy_nlb = [None, None, None]
    for n, lib in t9.items():
        if lib['_type'] == 'decay':
            decay_nlb.append(n)
        elif lib['_subtype'] == 'activation_products':
            xsfpy_nlb[0] = n
        elif lib['_subtype'] == 'actinides':
            xsfpy_nlb[1] = n
        elif lib['_subtype'] == 'fission_products':
            xsfpy_nlb[2] = n
    decay_nlb.sort()
    return tuple(decay_nlb), tuple(xsfpy_nlb)


def make_tape9(nucs, xscache=None, nlb=(201, 202, 203)):
    """Make a TAPE9 dict with data for a given list of nucs using data from
    a given data source.

    Parameters
    ----------
    nucs : iterable of ints, optional
        Set of nuclides in any format.

    Returns
    -------
    tape9: dict
        A full TAPE9 nested structure inside a dict. Keys 1, 2, and 3 correspond
        to decay decks. 219 is the activation products deck. 220 is the
        actinides deck. 221 is the fission product yield deck.
    """
    # build decay decks
    decay_file = StringIO(decay_tape9.decay_tape9)
    decay = parse_tape9(decay_file)

    if xscache is None:
        xscache = cache.XSCache()
    nucs = {nucname.id(nuc) for nuc in nucs}
    xsfpys = xslibs(nucs=nucs, xscache=xscache, nlb=nlb)
    tape9 = merge_tape9([decay, xsfpys])
    return tape9
