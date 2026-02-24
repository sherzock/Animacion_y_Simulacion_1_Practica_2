using MathNet.Numerics.Providers.LinearAlgebra;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using static UnityEditor.PlayerSettings;
using DenseMatrixXD = MathNet.Numerics.LinearAlgebra.Double.DenseMatrix;
using DenseVectorXD = MathNet.Numerics.LinearAlgebra.Double.DenseVector;
using MatrixXD = MathNet.Numerics.LinearAlgebra.Matrix<double>;
using VectorXD = MathNet.Numerics.LinearAlgebra.Vector<double>;

/// <summary>
/// Basic point constraint between two rigid bodies.
/// </summary>
public class PointConstraint : MonoBehaviour, IConstraint
{
    /// <summary>
    /// Default constructor. All zero. 
    /// </summary>
    public PointConstraint()
    {
        Manager = null;
    }

    #region EditorVariables

    public float Stiffness;

    public RigidBody bodyA;
    public RigidBody bodyB;

    #endregion

    #region OtherVariables

    int index;
    private PhysicsManager Manager;

    protected Vector3 pointA;
    protected Vector3 pointB;

    #endregion

    #region MonoBehaviour

    // Update is called once per frame
    void Update()
    {
        // Compute the average position
        Vector3 posA = (bodyA != null) ? bodyA.PointLocalToGlobal(pointA) : pointA;
        Vector3 posB = (bodyB != null) ? bodyB.PointLocalToGlobal(pointB) : pointB;
        Vector3 pos = 0.5f * (posA + posB);

        // Apply the position
        Transform xform = GetComponent<Transform>();
        xform.position = pos;
    }

    #endregion

    #region IConstraint

    public void Initialize(int ind, PhysicsManager m)
    {
        index = ind;
        Manager = m;

        // Initialize local positions. We assume that the object is connected to a Sphere mesh.
        Transform xform = GetComponent<Transform>();
        if (xform == null)
        {
            System.Console.WriteLine("[ERROR] Couldn't find any transform to the constraint");
        }
        else
        {
            System.Console.WriteLine("[TRACE] Succesfully found transform connected to the constraint");
        }

        // Initialize kinematics
        Vector3 pos = xform.position;

        // Local positions on objects
        pointA = (bodyA != null) ? bodyA.PointGlobalToLocal(pos) : pos;
        pointB = (bodyB != null) ? bodyB.PointGlobalToLocal(pos) : pos;

    }

    public int GetNumConstraints()
    {
        return 3;
    }

    public void GetConstraints(VectorXD c)
    {
        // TO BE COMPLETED
        Vector3 posA = (bodyA != null) ? bodyA.PointLocalToGlobal(pointA) : pointA;
        Vector3 posB = (bodyB != null) ? bodyB.PointLocalToGlobal(pointB) : pointB;

        float cx = (posA[0] - posB[0]);
        float cy = (posA[1] - posB[1]);
        float cz = (posA[2] - posB[2]);
         c = Utils.ToVectorXD(new Vector3(cx, cy, cz));
    }

    public void GetConstraintJacobian(MatrixXD dcdx)
    {
        // TO BE COMPLETED
    }

    public void GetForce(VectorXD force)
    {

        Vector3 posA = (bodyA != null) ? bodyA.PointLocalToGlobal(pointA) : pointA;
        Vector3 posB = (bodyB != null) ? bodyB.PointLocalToGlobal(pointB) : pointB;

        // TO BE COMPLETED
        float cx = (posA[0] - posB[0]);
        float cy = (posA[1] - posB[1]);
        float cz = (posA[2] - posB[2]);
        VectorXD c = Utils.ToVectorXD(new Vector3(cx, cy, cz));
        if (bodyA != null)
        {
            MatrixXD dcdxa = - Utils.Skew(posA - bodyA.m_pos);

            VectorXD fax = -Stiffness * dcdxa.Transpose() * c;
            VectorXD faTheta = -Stiffness * dcdxa.Transpose() * c;
        
            force.SetSubVector(bodyA.index, 3, force.SubVector(bodyA.index, 3) + fax);
            force.SetSubVector(bodyA.index + 3, 3, force.SubVector(bodyA.index + 3, 3) + faTheta);
        }

        if (bodyB != null)
        {
            MatrixXD dcdxb = -Utils.Skew(posB - bodyB.m_pos);

            VectorXD fbx = -Stiffness * dcdxb.Transpose() * c;
            VectorXD fbTheta = -Stiffness * dcdxb.Transpose() * c;

            force.SetSubVector(bodyB.index, 3, force.SubVector(bodyB.index, 3) + fbx);
            force.SetSubVector(bodyB.index + 3, 3, force.SubVector(bodyB.index + 3, 3) + fbTheta);
        }


    }

    public void GetForceJacobian(MatrixXD dFdx, MatrixXD dFdv)
    {
        // TO BE COMPLETED
    }

    #endregion

}
