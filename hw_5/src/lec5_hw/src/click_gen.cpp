#include "lec5_hw/visualizer.hpp"
#include "lec5_hw/trajectory.hpp"

#include <ros/ros.h>
#include <geometry_msgs/Point.h>
#include <geometry_msgs/PoseStamped.h>

#include <cmath>
#include <iostream>
#include <vector>

struct Config
{
    std::string targetTopic;
    double clickHeight;
    std::vector<double> initialVel;
    std::vector<double> initialAcc;
    std::vector<double> terminalVel;
    std::vector<double> terminalAcc;
    double allocationSpeed;
    double allocationAcc;
    int maxPieceNum;

    Config(const ros::NodeHandle &nh_priv)
    {
        nh_priv.getParam("TargetTopic", targetTopic);
        nh_priv.getParam("ClickHeight", clickHeight);
        nh_priv.getParam("InitialVel", initialVel);
        nh_priv.getParam("InitialAcc", initialAcc);
        nh_priv.getParam("TerminalVel", terminalVel);
        nh_priv.getParam("TerminalAcc", terminalAcc);
        nh_priv.getParam("AllocationSpeed", allocationSpeed);
        nh_priv.getParam("AllocationAcc", allocationAcc);
        nh_priv.getParam("MaxPieceNum", maxPieceNum);
    }
};

double timeTrapzVel(const double dist,
                    const double vel,
                    const double acc)
{
    const double t = vel / acc;
    const double d = 0.5 * acc * t * t;

    if (dist < d + d)
    {
        return 2.0 * sqrt(dist / acc);
    }
    else
    {
        return 2.0 * t + (dist - 2.0 * d) / vel;
    }
}
/*
//get M
Eigen::MatrixXd get_M(
    const Eigen::VectorXd &timeAllocationVector,
    const int pieceNum)
{
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(6*pieceNum,6*pieceNum);
    
    for(int i=0;i<pieceNum;i++)
    {
        Eigen::MatrixXd M_k = Eigen::MatrixXd::Zero(6,6);
        M_k(0,0)=1;
        M_k(1,1)=1;
        M_k(2,2)=2;
        for(int j=0;j<6;j++)
        {
            M_k(3,j)=pow(timeAllocationVector(i),j);
        }
        for(int j=0;j<6;j++)
        {
            M_k(4,j)=j*pow(timeAllocationVector(i),(j-1));
        }
        M_k(5,5)=20*pow(timeAllocationVector(i),3);
        M_k(5,4)=12*pow(timeAllocationVector(i),2);
        M_k(5,3)=6*pow(timeAllocationVector(i),1);
        M_k(5,2)=2*pow(timeAllocationVector(i),0);

        M.block(6*i,6*i,6,6)= M_k;
    }

    return M;
    
}

//get Q
Eigen::MatrixXd get_Q(
    const Eigen::VectorXd &timeAllocationVector,
    const int pieceNum)
{
    Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(6*pieceNum,6*pieceNum);
    
    for(int i=0;i<pieceNum;i++)
    {
        Eigen::MatrixXd Q_k = Eigen::MatrixXd::Zero(3,3);
        Q_k(3,3)=18*pow(timeAllocationVector(i),2);
        Q_k(3,4)=24*pow(timeAllocationVector(i),3);
        Q_k(3,5)=30*pow(timeAllocationVector(i),4);
        Q_k(4,3)=24*pow(timeAllocationVector(i),3);
        Q_k(4,4)=48*pow(timeAllocationVector(i),4);
        Q_k(4,5)=72*pow(timeAllocationVector(i),5);
        Q_k(5,3)=30*pow(timeAllocationVector(i),4);
        Q_k(5,4)=72*pow(timeAllocationVector(i),5);
        Q_k(5,5)=120*pow(timeAllocationVector(i),6);

        Q.block(3*(2*i+1),3*(2*i+1),3,3)= Q_k;
    }

    return Q;
    
}

//get Ct
Eigen::MatrixXd get_Ct(const int pieceNum)
{
    int ct_row=6*pieceNum;
    int ct_col=3*pieceNum+3;
    Eigen::MatrixXd Ct=Eigen::MatrixXd::Zero(ct_row,ct_col);
    for(int i=0;i<3;i++)
    {
        Ct(i,i)=1;
    }
    for(int i=0;i<3;i++)
    {
        Ct(6*pieceNum-3+i,2+pieceNum+i)=1;
    }
    if(pieceNum>1)
    {
        for(int i=1;i<pieceNum;i++)
        {
            Ct(6*i+1-4,3+i-1)=1;
            Ct(6*i+2-4,pieceNum+5+2*(i-1)+0)=1;
            Ct(6*i+3-4,pieceNum+5+2*(i-1)+1)=1;

            Ct(6*i+4-4,3+i-1)=1;
            Ct(6*i+5-4,pieceNum+5+2*(i-1)+0)=1;
            Ct(6*i+6-4,pieceNum+5+2*(i-1)+1)=1;
        }
    }
    return Ct;

}


*/

void minimumJerkTrajGen(
    // Inputs:
    const int pieceNum,
    const Eigen::Vector3d &initialPos,
    const Eigen::Vector3d &initialVel,
    const Eigen::Vector3d &initialAcc,
    const Eigen::Vector3d &terminalPos,
    const Eigen::Vector3d &terminalVel,
    const Eigen::Vector3d &terminalAcc,
    const Eigen::Matrix3Xd &intermediatePositions,
    const Eigen::VectorXd &timeAllocationVector,
    // Outputs:
    Eigen::MatrixX3d &coefficientMatrix)
{
    // coefficientMatrix is a matrix with 6*piece num rows and 3 columes
    // As for a polynomial c0+c1*t+c2*t^2+c3*t^3+c4*t^4+c5*t^5,
    // each 6*3 sub-block of coefficientMatrix is
    // --              --
    // | c0_x c0_y c0_z |
    // | c1_x c1_y c1_z |
    // | c2_x c2_y c2_z |
    // | c3_x c3_y c3_z |
    // | c4_x c4_y c4_z |
    // | c5_x c5_y c5_z |
    // --              --
    // Please computed coefficientMatrix of the minimum-jerk trajectory
    // in this function

    // ------------------------ Put your solution below ------------------------

    /*
    Eigen::MatrixXd Q,M,Ct,C,R,inv_M,R_pp,R_fp;
    Q=get_Q(timeAllocationVector,pieceNum);
    M=get_M(timeAllocationVector,pieceNum);
    Ct=get_Ct(pieceNum);

    C=Ct.transpose();
    inv_M=M.inverse();
    R=C*(inv_M.transpose())*Q*inv_M*Ct;

    R_pp=R.block(pieceNum+5,pieceNum+5,2*(pieceNum-1),2*(pieceNum-1));
    R_fp=R.block(0,pieceNum+5,pieceNum+5,2*(pieceNum-1));

    Eigen::MatrixXd dFx=Eigen::MatrixXd::Zero(pieceNum+5,1);
    Eigen::MatrixXd dFy=Eigen::MatrixXd::Zero(pieceNum+5,1);
    Eigen::MatrixXd dFz=Eigen::MatrixXd::Zero(pieceNum+5,1);


    Eigen::MatrixXd dPx=-R_pp.inverse()*R_fp.transpose()*dFx;
    */
    
    Eigen::MatrixXd M=Eigen::MatrixXd::Zero(6*pieceNum,6*pieceNum);
    Eigen::MatrixXd b=Eigen::MatrixXd::Zero(6*pieceNum,3);
 
    Eigen::MatrixXd F_0(3,6);
    F_0 <<  1, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0,
            0, 0, 2, 0, 0, 0;
    M.block(0, 0, 3, 6) = F_0;
    b.block(0, 0, 3, 3) <<  initialPos(0), initialPos(1), initialPos(2), 
                            initialVel(0), initialVel(1), initialVel(2), 
                            initialAcc(0), initialAcc(1), initialAcc(2);
    
    for (int i=1; i < pieceNum ; i++){
        double t(timeAllocationVector(i-1));
        Eigen::MatrixXd F_i(6,6), E_i(6,6);
        Eigen::Vector3d D_i(intermediatePositions.transpose().row(i-1));
        E_i <<  1, t, pow(t, 2), pow(t, 3), pow(t, 4), pow(t, 5),
                1, t, pow(t, 2), pow(t, 3), pow(t, 4), pow(t, 5),
                0, 1, 2*t, 3*pow(t, 2), 4*pow(t, 3), 5*pow(t, 4),
                0, 0, 2, 6*t, 12*pow(t, 2), 20*pow(t, 3),
                0, 0, 0, 6, 24*t, 60*pow(t, 2),
                0, 0, 0, 0, 24, 120*t;
        F_i <<  0, 0, 0, 0, 0, 0,
                -1, 0, 0, 0, 0, 0,
                0, -1, 0, 0, 0, 0,
                0, 0, -2, 0, 0, 0,
                0, 0, 0, -6, 0, 0,
                0, 0, 0, 0, -24, 0;
        int j = 6 * (i-1);
        M.block(3+6*(i-1), j + 6, 6, 6) = F_i;
        M.block(3+6*(i-1), j, 6, 6) = E_i;
        b.block(3+6*(i-1), 0, 6, 3) <<  D_i(0), D_i(1), D_i(2),
                                        0, 0, 0,
                                        0, 0, 0,
                                        0, 0, 0,
                                        0, 0, 0,
                                        0, 0, 0;
    }
    double t(timeAllocationVector(pieceNum-1));
    Eigen::MatrixXd E_M(3,6);
    E_M <<  1, t, pow(t, 2), pow(t, 3), pow(t, 4), pow(t, 5),
                0, 1, 2*t, 3*pow(t, 2), 4*pow(t, 3), 5*pow(t, 4),
                0, 0, 2, 6*t, 12*pow(t, 2), 20*pow(t, 3);
    M.block(6 * pieceNum -3, 6 * (pieceNum - 1), 3, 6) = E_M;
    b.block(6 * pieceNum -3, 0, 3, 3) <<    terminalPos(0), terminalPos(1), terminalPos(2), 
                                            terminalVel(0), terminalVel(1), terminalVel(2), 
                                            terminalAcc(0), terminalAcc(1), terminalAcc(2);
                                        
 
    coefficientMatrix = M.inverse() * b;
    
    // ------------------------ Put your solution above ------------------------
}

class ClickGen
{
private:
    Config config;

    ros::NodeHandle nh;
    ros::Subscriber targetSub;

    Visualizer visualizer;

    Eigen::Matrix3Xd positions;
    Eigen::VectorXd times;
    int positionNum;
    Trajectory<5> traj;

public:
    ClickGen(const Config &conf,
             ros::NodeHandle &nh_)
        : config(conf),
          nh(nh_),
          visualizer(nh),
          positions(3, config.maxPieceNum + 1),
          times(config.maxPieceNum),
          positionNum(0)
    {
        targetSub = nh.subscribe(config.targetTopic, 1,
                                 &ClickGen::targetCallBack, this,
                                 ros::TransportHints().tcpNoDelay());
    }

    void targetCallBack(const geometry_msgs::PoseStamped::ConstPtr &msg)
    {
        if (positionNum > config.maxPieceNum)
        {
            positionNum = 0;
            traj.clear();
        }

        positions(0, positionNum) = msg->pose.position.x;
        positions(1, positionNum) = msg->pose.position.y;
        positions(2, positionNum) = std::fabs(msg->pose.orientation.z) * config.clickHeight;

        if (positionNum > 0)
        {
            const double dist = (positions.col(positionNum) - positions.col(positionNum - 1)).norm();
            times(positionNum - 1) = timeTrapzVel(dist, config.allocationSpeed, config.allocationAcc);
        }

        ++positionNum;

        if (positionNum > 1)
        {
            const int pieceNum = positionNum - 1;
            const Eigen::Vector3d initialPos = positions.col(0);
            const Eigen::Vector3d initialVel(config.initialVel[0], config.initialVel[1], config.initialVel[2]);
            const Eigen::Vector3d initialAcc(config.initialAcc[0], config.initialAcc[1], config.initialAcc[2]);
            const Eigen::Vector3d terminalPos = positions.col(pieceNum);
            const Eigen::Vector3d terminalVel(config.terminalVel[0], config.terminalVel[1], config.terminalVel[2]);
            const Eigen::Vector3d terminalAcc(config.terminalAcc[0], config.terminalAcc[1], config.terminalAcc[2]);
            const Eigen::Matrix3Xd intermediatePositions = positions.middleCols(1, pieceNum - 1);
            const Eigen::VectorXd timeAllocationVector = times.head(pieceNum);

            Eigen::MatrixX3d coefficientMatrix = Eigen::MatrixXd::Zero(6 * pieceNum, 3);

            minimumJerkTrajGen(pieceNum,
                               initialPos, initialVel, initialAcc,
                               terminalPos, terminalVel, terminalAcc,
                               intermediatePositions,
                               timeAllocationVector,
                               coefficientMatrix);

            traj.clear();
            traj.reserve(pieceNum);
            for (int i = 0; i < pieceNum; i++)
            {
                traj.emplace_back(timeAllocationVector(i),
                                  coefficientMatrix.block<6, 3>(6 * i, 0).transpose().rowwise().reverse());
            }
        }

        visualizer.visualize(traj, positions.leftCols(positionNum));

        return;
    }
};

int main(int argc, char **argv)
{
    ros::init(argc, argv, "click_gen_node");
    ros::NodeHandle nh_;
    ClickGen clickGen(Config(ros::NodeHandle("~")), nh_);
    ros::spin();
    return 0;
}
