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

        VectorXD cv = Utils.ToVectorXD(posA - posB);
        c.SetSubVector(index, 3, cv);
    }

    private static void AddBlock(MatrixXD A, int r0, int c0, MatrixXD block3x3)
    {
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++)
                A[r0 + r, c0 + c] += block3x3[r, c];
    }
    public void GetConstraintJacobian(MatrixXD dcdx)
    {
        // TO BE COMPLETED
        int row = index;
        Vector3 posA = (bodyA != null) ? bodyA.PointLocalToGlobal(pointA) : pointA;
        Vector3 posB = (bodyB != null) ? bodyB.PointLocalToGlobal(pointB) : pointB;
        MatrixXD I3 = DenseMatrixXD.CreateIdentity(3);

        if (bodyA != null)
        {
            int colA = bodyA.index;
            AddBlock(dcdx, row, colA, I3);
            Vector3 rA = posA - bodyA.m_pos;
            AddBlock(dcdx, row, colA + 3, -Utils.Skew(rA));
        }

        if (bodyB != null)
        {
            int colB = bodyB.index;
            AddBlock(dcdx, row, colB, -I3);
            Vector3 rB = posB - bodyB.m_pos;
            AddBlock(dcdx, row, colB + 3, Utils.Skew(rB));
        }
    }

    public void GetForce(VectorXD force)
    {

        Vector3 posA = (bodyA != null) ? bodyA.PointLocalToGlobal(pointA) : pointA;
        Vector3 posB = (bodyB != null) ? bodyB.PointLocalToGlobal(pointB) : pointB;

        // TO BE COMPLETED
        Vector3 c3 = posA - posB;
        VectorXD c = Utils.ToVectorXD(c3);

        if (bodyA != null)
        {
            MatrixXD dcdxa = - Utils.Skew(posA - bodyA.m_pos);

            VectorXD fax = -Stiffness *  c;
            VectorXD faTheta = -Stiffness * (dcdxa.Transpose() * c);
        
            force.SetSubVector(bodyA.index, 3, force.SubVector(bodyA.index, 3) + fax);
            force.SetSubVector(bodyA.index + 3, 3, force.SubVector(bodyA.index + 3, 3) + faTheta);
        }

        if (bodyB != null)
        {
            MatrixXD dcdxb = Utils.Skew(posB - bodyB.m_pos);

            VectorXD fbx = Stiffness *  c;
            VectorXD fbTheta = -Stiffness * (dcdxb.Transpose() * c);

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
