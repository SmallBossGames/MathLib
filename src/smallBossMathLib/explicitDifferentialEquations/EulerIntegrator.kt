package smallBossMathLib.explicitDifferentialEquations

import smallBossMathLib.shared.IntegratorBase
import smallBossMathLib.shared.zeroSafetyNorm
import kotlin.math.sqrt


class EulerIntegrator(val startEvaluations: Int,
                      val accuracy: Double,
                      val minStep: Double,
                      val maxStep: Double) : IntegratorBase() {
    fun integrate(
        t0: Double,
        y0: DoubleArray,
        t: Double,
        outY: DoubleArray,
        equations: (y: DoubleArray, f: DoubleArray) -> Unit
    ) {
        executeStepHandlers(t0, y0)

        y0.copyInto(outY)

        var fCurrentBuffer = DoubleArray(y0.size)
        var fNextBuffer = DoubleArray(y0.size)
        var currentEvaluationsCount = 0

        val yNextBuffer = DoubleArray(y0.size)

        var time = t0
        var step = t / startEvaluations

        val vectorBuffer = DoubleArray(y0.size)

        val endTime = time+t

        equations(outY, fCurrentBuffer)

        while (time < endTime){
            for (i in yNextBuffer.indices){
                yNextBuffer[i] = outY[i] + step*fCurrentBuffer[i]
            }

            equations(yNextBuffer, fNextBuffer)

            for (i in vectorBuffer.indices){
                vectorBuffer[i] = fNextBuffer[i] - fCurrentBuffer[i]
            }

            val errNorm = 0.5 * step * zeroSafetyNorm(vectorBuffer, outY, 1e-7)
            val q = sqrt(accuracy / errNorm)

            if(q < 1){
                step = q * step / 1.1
                continue
            }
            else {
                time += step
                step = q * step / 1.1

                for (i in outY.indices){
                    outY[i] = yNextBuffer[i]
                }

                val temp = fCurrentBuffer
                fCurrentBuffer = fNextBuffer
                fNextBuffer = temp
            }
        }
    }

}