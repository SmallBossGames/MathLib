package smallBossMathLib.explicitDifferentialEquations

import smallBossMathLib.shared.Integrator
import smallBossMathLib.shared.StepInfo
import kotlin.math.abs

class RKM4Integrator(val evaluations: Int, val accuracy: Double) : Integrator() {
    fun integrate(
        t0: Double,
        y0: DoubleArray,
        t: Double,
        rVector:DoubleArray,
        outY: DoubleArray,
        equations: (y: DoubleArray, outF: DoubleArray) -> Unit
    ) {
        y0.copyInto(outY)

        val fCurrentBuffer = DoubleArray(y0.size)
        val yNextBuffer = DoubleArray(y0.size)
        val vectorBuffer1 = DoubleArray(y0.size)

        val k = Array(5){ DoubleArray(y0.size) }

        var time = t0
        var step = t / evaluations

        var isLowStepSizeReached = isLowStepSizeReached(step)
        var isHighStepSizeReached = isHighStepSizeReached(step)

        val endTime = t0 + t

        val stepInfo = StepInfo(t0, y0, isLowStepSizeReached, isHighStepSizeReached)

        executeStepHandlers(stepInfo)

        mainLoop@ while (time < endTime) {
            step = normalizeStep(step, time, endTime)

            equations(outY, fCurrentBuffer)

            for (i in fCurrentBuffer.indices) {
                k[0][i] = step * fCurrentBuffer[i]
                yNextBuffer[i] = outY[i] + 1.0 / 3.0 * k[0][i]
            }

            equations(yNextBuffer, fCurrentBuffer)

            for (i in fCurrentBuffer.indices) {
                k[1][i] = step * fCurrentBuffer[i]
                yNextBuffer[i] = outY[i] + 1.0 / 6.0 * k[0][i] + 1.0 / 6.0 * k[1][i]
            }

            equations(yNextBuffer, fCurrentBuffer)

            for (i in fCurrentBuffer.indices) {
                k[2][i] = step * fCurrentBuffer[i]
                yNextBuffer[i] = outY[i] + 1.0 / 8.0 * k[0][i] + 3.0 / 8.0 * k[2][i]
            }

            equations(yNextBuffer, fCurrentBuffer)

            for (i in fCurrentBuffer.indices) {
                k[3][i] = step * fCurrentBuffer[i]
                yNextBuffer[i] = outY[i] + 1.0 / 2.0 * k[0][i] - 3.0 / 2.0 * k[2][i] + 2.0 * k[3][i]
            }

            equations(yNextBuffer, fCurrentBuffer)

            for (i in fCurrentBuffer.indices) {
                k[4][i] = step * fCurrentBuffer[i]
            }

            for (i in vectorBuffer1.indices){
                vectorBuffer1[i] = abs((2.0*k[0][i] - 9.0*k[2][i] + 8.0*k[3][i] - k[4][i]) / 30.0)
            }

            if(!isLowStepSizeReached){
                for(i in vectorBuffer1.indices){
                    if(abs(outY[i]) > rVector[i] && vectorBuffer1[i] >= accuracy*abs(outY[i])){
                        step /= 2.0
                        isLowStepSizeReached = isLowStepSizeReached(step)
                        continue@mainLoop
                    }
                }
            }

            for (i in outY.indices){
                outY[i] += 1.0/6.0*k[0][i] + 2.0/3.0*k[3][i] + 1.0/6.0*k[4][i]
            }

            time += step

            stepInfo.set(time, outY, isLowStepSizeReached, isHighStepSizeReached)

            executeStepHandlers(stepInfo)

            if (!isHighStepSizeReached) {
                for (i in vectorBuffer1.indices) {
                    if(vectorBuffer1[i] > accuracy*abs(outY[i]) / 32.0){
                        continue@mainLoop
                    }
                }
                step *= 2.0
                isHighStepSizeReached = isHighStepSizeReached(step)
            }
        }
    }
}