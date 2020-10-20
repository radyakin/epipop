
program define epipop_siz_at
	version 16.0
	syntax , [clear]
	
	tempname M
	matrix `M'= ///
	0,1.580115e-06,3.74938e-07 \ ///
	1,1.580115e-06,3.74938e-07 \ ///
	2,1.580115e-06,3.74938e-07 \ ///
	3,1.580115e-06,3.74938e-07 \ ///
	4,1.580115e-06,3.74938e-07 \ ///
	5,1.580115e-06,3.74938e-07 \ ///
	6,1.580115e-06,3.74938e-07 \ ///
	7,1.580115e-06,3.74938e-07 \ ///
	8,1.580115e-06,3.74938e-07 \ ///
	9,1.580115e-06,3.74938e-07 \ ///
	10,7.81733e-07,7.31053e-07 \ ///
	11,7.81733e-07,7.31053e-07 \ ///
	12,7.81733e-07,7.31053e-07 \ ///
	13,7.81733e-07,7.31053e-07 \ ///
	14,7.81733e-07,7.31053e-07 \ ///
	15,7.81733e-07,7.31053e-07 \ ///
	16,7.81733e-07,7.31053e-07 \ ///
	17,7.81733e-07,7.31053e-07 \ ///
	18,7.81733e-07,7.31053e-07 \ ///
	19,7.81733e-07,7.31053e-07 \ ///
	20,8.953379e-06,7.834254e-06 \ ///
	21,8.953379e-06,7.834254e-06 \ ///
	22,8.953379e-06,7.834254e-06 \ ///
	23,8.953379e-06,7.834254e-06 \ ///
	24,8.953379e-06,7.834254e-06 \ ///
	25,8.953379e-06,7.834254e-06 \ ///
	26,8.953379e-06,7.834254e-06 \ ///
	27,8.953379e-06,7.834254e-06 \ ///
	28,8.953379e-06,7.834254e-06 \ ///
	29,8.953379e-06,7.834254e-06 \ ///
	30,.00005449847,.00002406607 \ ///
	31,.00005449847,.00002406607 \ ///
	32,.00005449847,.00002406607 \ ///
	33,.00005449847,.00002406607 \ ///
	34,.00005449847,.00002406607 \ ///
	35,.00005449847,.00002406607 \ ///
	36,.00005449847,.00002406607 \ ///
	37,.00005449847,.00002406607 \ ///
	38,.00005449847,.00002406607 \ ///
	39,.00005449847,.00002406607 \ ///
	40,.0001718814,.00009954241 \ ///
	41,.0001718814,.00009954241 \ ///
	42,.0001718814,.00009954241 \ ///
	43,.0001718814,.00009954241 \ ///
	44,.0001718814,.00009954241 \ ///
	45,.0001718814,.00009954241 \ ///
	46,.0001718814,.00009954241 \ ///
	47,.0001718814,.00009954241 \ ///
	48,.0001718814,.00009954241 \ ///
	49,.0001718814,.00009954241 \ ///
	50,.0004394879,.0002251994 \ ///
	51,.0004394879,.0002251994 \ ///
	52,.0004394879,.0002251994 \ ///
	53,.0004394879,.0002251994 \ ///
	54,.0004394879,.0002251994 \ ///
	55,.0004394879,.0002251994 \ ///
	56,.0004394879,.0002251994 \ ///
	57,.0004394879,.0002251994 \ ///
	58,.0004394879,.0002251994 \ ///
	59,.0004394879,.0002251994 \ ///
	60,.00013349083,.00009159046 \ ///
	61,.00013349083,.00009159046 \ ///
	62,.00013349083,.00009159046 \ ///
	63,.00013349083,.00009159046 \ ///
	64,.00013349083,.00009159046 \ ///
	65,.00013349083,.00009159046 \ ///
	66,.00013349083,.00009159046 \ ///
	67,.00013349083,.00009159046 \ ///
	68,.00013349083,.00009159046 \ ///
	69,.00013349083,.00009159046 \ ///
	70,.00013349083,.00009159046 \ ///
	71,.00013349083,.00009159046 \ ///
	72,.00013349083,.00009159046 \ ///
	73,.00013349083,.00009159046 \ ///
	74,.00013349083,.00009159046 \ ///
	75,.00013349083,.00009159046 \ ///
	76,.00013349083,.00009159046 \ ///
	77,.00013349083,.00009159046 \ ///
	78,.00013349083,.00009159046 \ ///
	79,.00013349083,.00009159046 \ ///
	80,.00013349083,.00009159046 \ ///
	81,.00013349083,.00009159046 \ ///
	82,.00013349083,.00009159046 \ ///
	83,.00013349083,.00009159046 \ ///
	84,.00013349083,.00009159046 \ ///
	85,.00013349083,.00009159046 \ ///
	86,.00013349083,.00009159046 \ ///
	87,.00013349083,.00009159046 \ ///
	88,.00013349083,.00009159046 \ ///
	89,.00013349083,.00009159046 \ ///
	90,.00013349083,.00009159046 \ ///
	91,.00013349083,.00009159046 \ ///
	92,.00013349083,.00009159046 \ ///
	93,.00013349083,.00009159046 \ ///
	94,.00013349083,.00009159046 \ ///
	95,.00013349083,.00009159046 \ ///
	96,.00013349083,.00009159046 \ ///
	97,.00013349083,.00009159046 \ ///
	98,.00013349083,.00009159046 \ ///
	99,.00013349083,.00009159046 \ ///
	100,.00013349083,.00009159046 \ ///
	101,.00013349083,.00009159046 \ ///
	102,.00013349083,.00009159046 \ ///
	103,.00013349083,.00009159046 \ ///
	104,.00013349083,.00009159046 \ ///
	105,.00013349083,.00009159046 \ ///
	106,.00013349083,.00009159046 \ ///
	107,.00013349083,.00009159046 \ ///
	108,.00013349083,.00009159046 \ ///
	109,.00013349083,.00009159046

	matrix colnames `M' = age male female
	`clear'
	svmat `M', names(col)
end

// end of file