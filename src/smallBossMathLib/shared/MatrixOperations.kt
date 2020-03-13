package smallBossMathLib.shared

import java.lang.Math.pow
import kotlin.math.abs
import kotlin.math.max
import kotlin.math.pow

const val rMin = 1e-14

fun findJacobiMatrix(inY:DoubleArray,
                     t: Double,
                     equations: (t: Double, inY: DoubleArray, outF: DoubleArray) -> Unit,
                     outJacobiMatrix: Array<DoubleArray>) {
    val tempY = inY.clone()
    val tempF = DoubleArray(tempY.size)
    val tempFWithDelta = DoubleArray(tempY.size)

    for (i in tempY.indices){

        val r = max(rMin, (rMin * tempY[i]))
        tempY[i] += r

        equations(t, inY, tempF)
        equations(t, tempY, tempFWithDelta)

        for (j in tempY.indices){
            outJacobiMatrix[j][i] = (tempFWithDelta[j] - tempF[j]) / r
        }

        tempY[i] = inY[i]
    }
}

fun multipleMatrix2D(constant: Double, matrix: Array<DoubleArray>, outMatrix: Array<DoubleArray>){
    for (i in outMatrix.indices)
        for (j in outMatrix[i].indices)
            outMatrix[i][j] = matrix[i][j] * constant
}

fun sumMatrix2D(m1: Array<DoubleArray>, m2: Array<DoubleArray>, outMatrix: Array<DoubleArray>){
    for (i in outMatrix.indices)
        for (j in outMatrix[i].indices)
            outMatrix[i][j] = m1[i][j] + m2[i][j]
}

fun subtractMatrix2D(m1: Array<DoubleArray>, m2: Array<DoubleArray>, outMatrix: Array<DoubleArray>){
    for (i in outMatrix.indices)
        for (j in outMatrix[i].indices)
            outMatrix[i][j] = m1[i][j] - m2[i][j]
}

fun createUnitMatrix(size: Int) : Array<DoubleArray>
{
    val e = Array(size){DoubleArray(size)}

    for (i in e.indices)
        e[i][i] = 1.0

    return e
}