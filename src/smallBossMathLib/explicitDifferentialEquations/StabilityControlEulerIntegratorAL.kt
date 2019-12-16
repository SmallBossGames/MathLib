package smallBossMathLib.explicitDifferentialEquations

import org.apache.commons.math3.analysis.solvers.UnivariateSolver
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations
import org.apache.commons.math3.ode.FirstOrderIntegrator
import org.apache.commons.math3.ode.events.EventHandler
import org.apache.commons.math3.ode.sampling.StepHandler

class StabilityControlEulerIntegratorAL : FirstOrderIntegrator {

    private var maxEvaluations: Int = 0
    private var evaluations: Int = 0

    override fun setMaxEvaluations(maxEvaluations: Int) {
        this.maxEvaluations = maxEvaluations
    }

    override fun getMaxEvaluations(): Int = maxEvaluations

    override fun getName(): String = "Stability control euler integrator"

    override fun getEvaluations(): Int = evaluations

    override fun getCurrentStepStart(): Double {
        TODO("not implemented") //To change body of created functions use File | Settings | File Templates.
    }

    override fun clearEventHandlers() {
        TODO("not implemented") //To change body of created functions use File | Settings | File Templates.
    }

    override fun getStepHandlers(): MutableCollection<StepHandler> {
        TODO("not implemented") //To change body of created functions use File | Settings | File Templates.
    }

    override fun addStepHandler(handler: StepHandler?) {
        TODO("not implemented") //To change body of created functions use File | Settings | File Templates.
    }

    override fun integrate(
        equations: FirstOrderDifferentialEquations?,
        t0: Double,
        y0: DoubleArray?,
        t: Double,
        y: DoubleArray?
    ): Double {
        TODO("not implemented") //To change body of created functions use File | Settings | File Templates.
    }

    override fun getEventHandlers(): MutableCollection<EventHandler> {
        TODO("not implemented") //To change body of created functions use File | Settings | File Templates.
    }

    override fun addEventHandler(
        handler: EventHandler?,
        maxCheckInterval: Double,
        convergence: Double,
        maxIterationCount: Int
    ) {
        TODO("not implemented") //To change body of created functions use File | Settings | File Templates.
    }

    override fun addEventHandler(
        handler: EventHandler?,
        maxCheckInterval: Double,
        convergence: Double,
        maxIterationCount: Int,
        solver: UnivariateSolver?
    ) {
        TODO("not implemented") //To change body of created functions use File | Settings | File Templates.
    }

    override fun clearStepHandlers() {
        TODO("not implemented") //To change body of created functions use File | Settings | File Templates.
    }

    override fun getCurrentSignedStepsize(): Double {
        TODO("not implemented") //To change body of created functions use File | Settings | File Templates.
    }

}