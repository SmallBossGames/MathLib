package smallBossMathLib.shared

import kotlin.math.abs
import kotlin.math.max

fun zeroSafetyNorm(vector: DoubleArray, yVector: DoubleArray, r: Double) : Double{
    var value = 0.0
    for (i in vector.indices)
        value = max(abs(vector[i]) / (abs(yVector[i]) + r), value)
    return value
}

fun sum(v1: DoubleArray, v2: DoubleArray, outV: DoubleArray){
    for (i in outV.indices)
        outV[i] = v1[i] + v2[i]
}