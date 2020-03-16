package smallBossMathLib.implicitDifferentialEquations

import smallBossMathLib.shared.Matrix2D
import smallBossMathLib.shared.findJacobiMatrix
import smallBossMathLib.shared.subtractMatrix2D

//const val a = 0.43586652150846

const val a = 0.292893219
const val p1 = 1.25
const val p2 = 0.75
const val beta = 0.666666667
const val alpha = -1.33333333

class MK22Integrator (val evaluations: Int, val accuracy: Double) {

    inline fun integrate(t0: Double,
                  y0: DoubleArray,
                  t: Double,
                  equations: (t: Double, inY: DoubleArray, outF: DoubleArray) -> Unit,
                  outY: DoubleArray){

        if(y0.size != outY.size)
            throw IllegalArgumentException()

        y0.copyInto(outY)

        val endTime = t0 + t
        val unitMatrix = Matrix2D.createUnitMatrix2D(y0.size)

        val matrixD = Matrix2D(y0.size)
        val k1 = DoubleArray(y0.size)
        val k2 = DoubleArray(y0.size)
        val yBuffer = DoubleArray(y0.size)
        val fBuffer = DoubleArray(y0.size)

        var step = t / evaluations
        var time = t0

        while (time < endTime){
            findJacobiMatrix(outY, t0, equations, matrixD)
            matrixD *= step
            matrixD *= a
            subtractMatrix2D(unitMatrix, matrixD, matrixD)
            matrixD.makeLU()

            for (i in yBuffer.indices){
                yBuffer[i] = outY[i]
            }
            equations(time, yBuffer, fBuffer)
            for (i in fBuffer.indices){
                fBuffer[i] = step * fBuffer[i]
            }
            matrixD.solveLU(fBuffer, k1)

            for (i in yBuffer.indices){
                yBuffer[i] = outY[i] + beta * k1[i]
            }
            equations(time, yBuffer, fBuffer)
            for (i in fBuffer.indices){
                fBuffer[i] = fBuffer[i] * step + alpha * k1[i]
            }
            matrixD.solveLU(fBuffer, k2)

            for (i in yBuffer.indices){
                yBuffer[i] = outY[i] + p1*k1[i] + p2*k2[i]
            }

            for (i in outY.indices){
                outY[i] = yBuffer[i]
            }
        }
    }


}