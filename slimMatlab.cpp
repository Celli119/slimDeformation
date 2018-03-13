//
// Created by girinon on 22/10/17.
//

#include <igl/doublearea.h>
#include <igl/flipped_triangles.h>
#include <igl/EPS.h>
#include <igl/slice.h>
#include <igl/euler_characteristic.h>
#include <igl/C_STR.h>
#include <igl/slim.h>
#include <igl/is_edge_manifold.h>
#include <igl/harmonic.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/boundary_loop.h>
#include <igl/components.h>
#include <igl/matlab/MexStream.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/validate_arg.h>
#include <igl/lim/lim.h>

void mexFunction(
        int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    const auto mexErrMsgTxt = [](const bool v ,const char * msg)
    {
        igl::matlab::mexErrMsgTxt(v,msg);
    };
    igl::matlab::MexStream mout;
    std::streambuf *outbuf = std::cout.rdbuf(&mout);
    mexErrMsgTxt(nrhs >= 4,"Four arguments expected");

    Eigen::MatrixX3d V,U0,U,bc;
    Eigen::MatrixXi F;
    Eigen::VectorXi b;

    int iters = 100;
    double p = 1e5;
    std::string solver;
    enum solverType {SLIM, LIM};
    solverType type(LIM);

    igl::matlab::parse_rhs_double(prhs+0,V);
    igl::matlab::parse_rhs_index(prhs+1,F);
    igl::matlab::parse_rhs_index(prhs+2,b);
    igl::matlab::parse_rhs_double(prhs+3,bc);

//    typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> MatrixXb;
    Eigen::MatrixXd dimCons;
    dimCons.resize(bc.rows(), bc.cols());
    dimCons.fill(true);

    std::cout << "number of elements:"<< F.rows() << std::endl;
    std::cout << "number of vertices:"<< V.rows() << std::endl;
    int i = 4;
    while(i<nrhs)
    {
        mexErrMsgTxt(mxIsChar(prhs[i]),"Parameter names should be strings");
        // Cast to char
        const char * name = mxArrayToString(prhs[i]);
        if(strcmp("P",name) == 0)
        {
            igl::matlab::validate_arg_scalar(i,nrhs,prhs,name);
            igl::matlab::validate_arg_double(i,nrhs,prhs,name);
            p = (double)*mxGetPr(prhs[++i]);
        }else if(strcmp("Iters",name) == 0)
        {
            igl::matlab::validate_arg_scalar(i,nrhs,prhs,name);
            igl::matlab::validate_arg_double(i,nrhs,prhs,name);
            iters = (double)*mxGetPr(prhs[++i]);
        }else if(strcmp("Solver",name) == 0)
        {
            igl::matlab::validate_arg_char(i,nrhs,prhs,name);
            solver = mxArrayToString(prhs[++i]);
            if (std::strcmp(solver.data(), "LIM") == 0){
                type = LIM;
            }else if (std::strcmp(solver.data(), "SLIM") == 0){
                type = SLIM;
            }
            else{
                type = LIM;
            }

        }else if(strcmp("DimCons",name) == 0)
        {
            igl::matlab::validate_arg_double(i,nrhs,prhs,name);
            i++;
            igl::matlab::parse_rhs_double(prhs + i, dimCons);
        }else
        {
            mexErrMsgTxt(false,C_STR("Unknown parameter: "<<name));
        }
        i++;
    }


    auto V_0 = V;
    if (type == SLIM){
        igl::SLIMData sData;
        sData.exp_factor = 5.0;
        slim_precompute(V,F,V_0,sData,igl::SLIMData::EXP_CONFORMAL,b,bc,dimCons,p);
        slim_solve(sData,iters);
        U = sData.V_o;
    }
    else
    {
        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(bc.rows() * 3);
        std::vector<double> constraintTargetsTmp; constraintTargetsTmp.reserve(V_0.rows() * 3);
        int maxRow = -1;
        for (int iConstrain = 0; iConstrain < b.rows(); ++iConstrain)
        {
            auto iVertex = b(iConstrain);
            for (int iDim = 0; iDim < 3; ++iDim) {
                if (dimCons(iVertex, iDim)){
                    auto idConstrain = 3 * iVertex + iDim;
                    if(maxRow < idConstrain)
                        maxRow = idConstrain;
                    constraintTargetsTmp.push_back(bc(iConstrain, iDim));
                    triplets.emplace_back(idConstrain,idConstrain,1);
                }
            }
        }
        auto constraintTargets = Eigen::VectorXd::Map(constraintTargetsTmp.data(), constraintTargetsTmp.size());

        Eigen::SparseMatrix<double> matConstraints;
        matConstraints.resize(maxRow + 1, V_0.rows() * 3);
        matConstraints.setFromTriplets(triplets.begin(), triplets.end());
        auto energy = igl::lim::Green;
        std::cout << matConstraints.size() << std::endl;
        U = V_0;
        std::vector<int> borderVertices;
        Eigen::VectorXd gradients;
        auto state = ComputeLIM(
                U,
                V_0,
                F,
                borderVertices,
                gradients,
                matConstraints,
                constraintTargets,
                (EnergyType)energy,
                1e-10,
                iters,
                true,
                true,
                true,
                true,
                -1,
                -1
        );

    }



    switch(nlhs)
    {
        default:
        {
            mexErrMsgTxt(false,"Too many output parameters.");
        }
        case 2:
        {
            igl::matlab::prepare_lhs_double(U0,plhs+1);
            // Fall through
        }
        case 1:
        {
            igl::matlab::prepare_lhs_double(U,plhs+0);
            // Fall through
        }
        case 0: break;
    }
    // Restore the std stream buffer Important!
    std::cout.rdbuf(outbuf);

}


