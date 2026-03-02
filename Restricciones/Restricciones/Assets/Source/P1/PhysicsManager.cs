using UnityEngine;
using System.Collections;
using System.Collections.Generic;

using VectorXD = MathNet.Numerics.LinearAlgebra.Vector<double>;
using MatrixXD = MathNet.Numerics.LinearAlgebra.Matrix<double>;
using DenseVectorXD = MathNet.Numerics.LinearAlgebra.Double.DenseVector;
using DenseMatrixXD = MathNet.Numerics.LinearAlgebra.Double.DenseMatrix;

/// <summary>
/// Basic physics manager capable of simulating a given ISimulable
/// implementation using diverse integration methods: explicit,
/// implicit, Verlet and semi-implicit.
/// </summary>
public class PhysicsManager : MonoBehaviour 
{
	/// <summary>
	/// Default constructor. Zero all. 
	/// </summary>
	public PhysicsManager()
	{
		Paused = true;
		TimeStep = 0.01f;
		Gravity = new Vector3 (0.0f, -9.81f, 0.0f);
		IntegrationMethod = Integration.Symplectic;
	}

	/// <summary>
	/// Integration method.
	/// </summary>
	public enum Integration
	{
		Symplectic = 1,
        Implicit = 2,
        SymplecticConstraints = 3,
    };

	#region InEditorVariables

	public bool Paused;
	public float TimeStep;
    public Vector3 Gravity;
    public List<GameObject> SimObjects;
    public List<GameObject> Constraints;
    public Integration IntegrationMethod;

    #endregion

    #region OtherVariables

    private List<ISimulable> m_objs;
    private List<IConstraint> m_constraints;
    private int m_numDoFs;
    private int m_numConstraints;

    #endregion

    #region MonoBehaviour

    public void Start()
    {
        //Parse the simulable objects and initialize their state indices
        m_numDoFs = 0;
        m_objs = new List<ISimulable>(SimObjects.Count);

        foreach (GameObject obj in SimObjects)
        {
            ISimulable simobj = obj.GetComponent<ISimulable>();
            if (simobj != null)
            {
                m_objs.Add(simobj);

                // Initialize simulable object
                simobj.Initialize(m_numDoFs, this);

                // Retrieve pos and vel size
                m_numDoFs += simobj.GetNumDoFs();
            }
        }

        //Parse the constraints
        m_numConstraints = 0;
        m_constraints = new List<IConstraint>(Constraints.Count);

        foreach (GameObject obj in Constraints)
        {
            IConstraint constraint = obj.GetComponent<IConstraint>();
            if (constraint != null)
            {
                m_constraints.Add(constraint);

                // Initialize constraint
                constraint.Initialize(m_numConstraints, this);

                // Retrieve the number of constraints
                m_numConstraints += constraint.GetNumConstraints();
            }
        }

    }

    public void Update()
	{
		if (Input.GetKeyUp (KeyCode.P))
			this.Paused = !this.Paused;

    }

    public void FixedUpdate()
    {
        if (this.Paused)
            return; // Not simulating

        // Select integration method
        switch (this.IntegrationMethod)
        {
            case Integration.Symplectic: this.stepSymplectic(); break;
            case Integration.Implicit: this.stepImplicit(); break;
            case Integration.SymplecticConstraints: this.stepSymplecticConstraints(); break;
            default:
                throw new System.Exception("[ERROR] Should never happen!");
        }
    }

    #endregion

    /// <summary>
    /// Performs a simulation step using Symplectic integration.
    /// </summary>
    private void stepSymplectic()
	{
        // TO BE COMPLETED
        VectorXD x = new DenseVectorXD(m_numDoFs);
        VectorXD v = new DenseVectorXD(m_numDoFs);
        VectorXD f = new DenseVectorXD(m_numDoFs);
        f.Clear();
        MatrixXD Minv = new DenseMatrixXD(m_numDoFs);
        Minv.Clear();

        foreach (ISimulable obj in m_objs)
        {
            obj.GetPosition(x);
            obj.GetVelocity(v);
            obj.GetForce(f);
            obj.GetMassInverse(Minv);
        }

        foreach(IConstraint constraint in m_constraints)
        {
            constraint.GetForce(f);
        }

        foreach (ISimulable obj in m_objs)
        {
            obj.FixVector(f);
            obj.FixMatrix(Minv);
        }

        v += TimeStep * (Minv * f);

        x = TimeStep * v;

        foreach (ISimulable obj in m_objs)
        {
            obj.AdvanceIncrementalPosition(x);
            obj.SetVelocity(v);
        }
    }

    /// <summary>
    /// Performs a simulation step using Implicit integration.
    /// </summary>
    private void stepImplicit()
    {
        // TO BE COMPLETED
        double h = TimeStep;

        VectorXD x = new DenseVectorXD(m_numDoFs);
        VectorXD v = new DenseVectorXD(m_numDoFs);

        foreach (ISimulable obj in m_objs)
        {
            obj.GetPosition(x);
            obj.GetVelocity(v);
        }



        VectorXD f = new DenseVectorXD(m_numDoFs); f.Clear();
        foreach (ISimulable obj in m_objs)
            obj.GetForce(f);

        foreach (IConstraint constraint in m_constraints)
        {
            constraint.GetForce(f);
        }

        MatrixXD Minv = new DenseMatrixXD(m_numDoFs);
        Minv.Clear();
        foreach (ISimulable obj in m_objs)
            obj.GetMassInverse(Minv);

        foreach (ISimulable obj in m_objs)
        {
            obj.FixVector(f);
            obj.FixMatrix(Minv);
        }

        MatrixXD dFdx = new DenseMatrixXD(m_numDoFs, m_numDoFs); dFdx.Clear();
        MatrixXD dFdv = new DenseMatrixXD(m_numDoFs, m_numDoFs); dFdv.Clear();
        
        foreach (ISimulable obj in m_objs)
            obj.GetForceJacobian(dFdx, dFdv);

        foreach (IConstraint constraint in m_constraints)
            constraint.GetForceJacobian(dFdx, dFdv);

        foreach (ISimulable obj in m_objs)
        {
            obj.FixMatrix(dFdx);
            obj.FixMatrix(dFdv);
        }

        MatrixXD I = DenseMatrixXD.CreateIdentity(m_numDoFs);

        MatrixXD A = I - h * (Minv * dFdv) - (h * h) * (Minv * dFdx);

        VectorXD b = (I - h * (Minv * dFdv)) * v + h * (Minv * f);

        foreach (ISimulable obj in m_objs)
        {
            obj.FixMatrix(A);
            obj.FixVector(b);
        }

        VectorXD vNew = A.Solve(b);

        VectorXD dx = h * vNew;

        foreach (ISimulable obj in m_objs)
        {
            obj.AdvanceIncrementalPosition(dx);
            obj.SetVelocity(vNew);
        }
    }

    /// <summary>
    /// Performs a simulation step using Symplectic integration with constrained dynamics.
    /// The constraints are treated as implicit
    /// </summary>
    private void stepSymplecticConstraints()
    {
        // TO BE COMPLETED
        double h = TimeStep;
        VectorXD x = new DenseVectorXD(m_numDoFs);
        VectorXD v = new DenseVectorXD(m_numDoFs);

        foreach (ISimulable obj in m_objs)
        {
            obj.GetPosition(x);
            obj.GetVelocity(v);
        }

        VectorXD f = new DenseVectorXD(m_numDoFs);
        f.Clear();
        foreach (ISimulable obj in m_objs)
            obj.GetForce(f);


        MatrixXD Minv = new DenseMatrixXD(m_numDoFs, m_numDoFs);
        Minv.Clear();
        foreach (ISimulable obj in m_objs)
            obj.GetMassInverse(Minv);

        foreach (ISimulable obj in m_objs)
        {
            obj.FixVector(f);
            obj.FixMatrix(Minv);
        }

        VectorXD vStar = v + h * (Minv * f);

        foreach (ISimulable obj in m_objs)
            obj.FixVector(vStar);

        VectorXD c = new DenseVectorXD(m_numConstraints);
        c.Clear();
        foreach (IConstraint con in m_constraints)
            con.GetConstraints(c);

        MatrixXD J = new DenseMatrixXD(m_numConstraints, m_numDoFs);
        J.Clear();
        foreach (IConstraint con in m_constraints)
            con.GetConstraintJacobian(J);

        foreach (ISimulable obj in m_objs)
        {
            obj.FixVector(c);
            obj.FixMatrix(J);
        }


        MatrixXD A = J * Minv * J.Transpose();
        double beta = 0.1;
        VectorXD rhs = -(J * vStar + (beta / h) * c);
        VectorXD lambda = A.Solve(rhs);

        VectorXD dv = Minv * (J.Transpose() * lambda);

        VectorXD vNew = vStar + dv;

        foreach (ISimulable obj in m_objs)
            obj.FixVector(vNew);

        VectorXD dx = h * vNew;

        foreach (ISimulable obj in m_objs)
        {
            obj.AdvanceIncrementalPosition(dx);
            obj.SetVelocity(vNew);
        }
    }

}
