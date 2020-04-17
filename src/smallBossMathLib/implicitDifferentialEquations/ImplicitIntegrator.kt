package smallBossMathLib.implicitDifferentialEquations

import smallBossMathLib.shared.Integrator

abstract class ImplicitIntegrator : Integrator() {
    private val stepHandlers = HashSet<(time:Double, y:DoubleArray, info: IImplicitMethodStepInfo) -> Unit>()

    protected fun executeStepHandlers(time:Double, y:DoubleArray, info: IImplicitMethodStepInfo) {
        for (handler in stepHandlers)
            handler(time, y, info)
    }

    fun addStepHandler(handler: (time:Double, y:DoubleArray, info: IImplicitMethodStepInfo) -> Unit) {
        stepHandlers += handler
    }

    fun removeStepHandler(handler: (time:Double, y:DoubleArray, info: IImplicitMethodStepInfo) -> Unit) {
        stepHandlers -= handler
    }
}