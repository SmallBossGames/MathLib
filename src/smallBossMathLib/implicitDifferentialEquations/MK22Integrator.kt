package smallBossMathLib.implicitDifferentialEquations

import smallBossMathLib.shared.Matrix2D
import smallBossMathLib.shared.findJacobiMatrix
import smallBossMathLib.shared.multipleMatrix2D
import smallBossMathLib.shared.subtractMatrix2D
import kotlin.math.abs
import kotlin.math.min
import kotlin.math.sqrt

//const val a = 0.43586652150846

const val a = 0.292893219
const val p1 = 1.25
const val p2 = 0.75
const val beta = 0.666666667
const val alpha = -1.33333333
const val accuracyCoefficient = (6.0*a - 5.0)/(4.0 - 8.0*a)

const val v = 1.0

class MK22Integrator (val evaluations: Int,
                      val freezeJacobiSteps: Int,
                      val stepSizeCoefficient: Double,
                      val accuracy: Double) {

    inline fun integrate(
        t0: Double,
        y0: DoubleArray,
        t: Double,
        outY: DoubleArray,
        equations: (t: Double, inY: DoubleArray, outF: DoubleArray) -> Unit) {

        if(y0.size != outY.size)
            throw IllegalArgumentException()

        y0.copyInto(outY)

        val endTime = t0 + t
        val unitMatrix = Matrix2D.createUnitMatrix2D(y0.size)

        val jacobiMatrix = Matrix2D(y0.size)
        val dMatrix = Matrix2D(y0.size)
        val k1 = DoubleArray(y0.size)
        val k2 = DoubleArray(y0.size)

        val vectorBuffer1 = DoubleArray(y0.size)
        val vectorBuffer2 = DoubleArray(y0.size)
        val matrixBuffer = Matrix2D(y0.size)

        var step = t / evaluations
        var time = t0
        var freezeSteps = 0
        var isNeedFindJacobi = true

        while (time < endTime){
            if(freezeSteps == 0 && isNeedFindJacobi){
                findJacobiMatrix(outY, t0, equations, jacobiMatrix)
            } else{
                isNeedFindJacobi = true
            }

            multipleMatrix2D(step * a, jacobiMatrix, dMatrix)
            subtractMatrix2D(unitMatrix, dMatrix, dMatrix)

            dMatrix.makeLU()

            for (i in vectorBuffer1.indices){
                vectorBuffer1[i] = outY[i]
            }
            equations(time, vectorBuffer1, vectorBuffer2)
            for (i in vectorBuffer1.indices){
                vectorBuffer2[i] = step * vectorBuffer2[i]
            }
            dMatrix.solveLU(vectorBuffer2, k1)

            for (i in vectorBuffer1.indices){
                vectorBuffer1[i] = outY[i] + beta * k1[i]
            }
            equations(time, vectorBuffer1, vectorBuffer2)
            for (i in vectorBuffer2.indices){
                vectorBuffer2[i] = step * vectorBuffer2[i] + alpha * k1[i]
            }
            dMatrix.solveLU(vectorBuffer2, k2)

            var e1Norm = 0.0
            for (i in vectorBuffer1.indices){
                val temp = accuracyCoefficient * (k2[i] + 1.0/3.0*k1[i])
                vectorBuffer1[i] = temp
                val norm = abs(temp) / (abs(outY[i]) + v)
                if (norm > e1Norm)
                    e1Norm = norm
            }

            val q1 = sqrt(8.0*accuracy / e1Norm)

            var e2Norm = 0.0
            val q2: Double
            if(q1 < 1.0) {
                dMatrix.inverseLU(matrixBuffer)
                matrixBuffer.multiply(vectorBuffer1, vectorBuffer2)
                for (i in vectorBuffer1.indices){
                    val norm = abs(vectorBuffer2[i]) / (abs(outY[i]) + v)
                    if (norm > e2Norm){
                        e2Norm = norm
                    }
                }
                q2 = sqrt(8.0*accuracy / e2Norm)
            } else {
                e2Norm = e1Norm
                q2 = q1
            }

            if(q2 < 1.0){
                step *= q2
                if(freezeSteps==0){
                    isNeedFindJacobi = false
                } else{
                    freezeSteps = 0
                }
                continue
            }

            for (i in outY.indices){
                outY[i] = outY[i] + p1*k1[i] + p2*k2[i]
            }

            time += step
            freezeSteps++

            val stepNew = min(q1, q2) * step

            if(!(freezeSteps < freezeJacobiSteps && stepNew < stepSizeCoefficient*step && e1Norm <= e2Norm)){
                freezeSteps = 0
                step = stepNew
            }
        }
    }
}