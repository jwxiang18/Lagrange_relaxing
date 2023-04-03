import gurobi.*;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

public class LR {
    Instance instance;
    int iterLimit = 100;
    int sameLimit = 3;
    double[] xlambda  ;
    double LB , UB = 0;
    List<Double> UBList = new LinkedList<>();
    List<Double> LBList = new LinkedList<>();
    List<Double> thetaList = new LinkedList<>();
    List<Double> stepList = new LinkedList<>();
    List<Double> gapList = new LinkedList<>();
    GRBModel model , modelRelax;
    GRBVar[][] vship;
    GRBVar[] vbuild;
    GRBConstr[] crleax;
    GRBConstr[] csupply;
    GRBConstr cbuildlimit;

    public LR(Instance instance){
        this.instance = instance;
    }

    public void main() throws GRBException {
        int same = 0;
        double theta = 1.0;
        double squareSum = 0.0;
        double step = 0.0;

        double sumDemand = Arrays.stream(instance.demand).sum();
        double selected_facility_supply = 0; // 实际的总运输量

        // 初始化模型
        init();

        // 初始化下界
        LB = relaxLB();

        xlambda = new double[instance.ncities];
        Arrays.fill(xlambda, 0.0);

        // 初始化上界
        for (int i = 0; i < instance.ncities; i++) {
            UB += Arrays.stream(instance.shipCost[i]).max().getAsDouble();
        }

        // 获得sum(cij * xij) 的表达式
        GRBLinExpr obj_shipCost = get_obj_shipCost();

        boolean isLRModel = false;

        for (int i = 0; i < iterLimit; i++) {
            if(!isLRModel){
                isLRModel = true;
                // 删除模型中松弛的约束
                for (GRBConstr grbConstr : crleax) {
                    model.remove(grbConstr);
                }
                crleax = new GRBConstr[instance.ncities];
            }

            GRBLinExpr obj_lagrangian = new GRBLinExpr();

            for (int j = 0; j < instance.ncities; j++) {
                for (int k = 0; k < instance.ncities; k++) {
                    obj_lagrangian.addTerm(-1*xlambda[j] , vship[k][j]);
                }
            }

            // 注意这里是一个常数项，在目标函数里没有添加，所以最后要单独添加
            double ud = 0;
            for (int j = 0; j < instance.ncities; j++) {
                ud += xlambda[j] * instance.demand[j];
            }
            obj_lagrangian.add(obj_shipCost);

            model.setObjective(obj_lagrangian , GRB.MINIMIZE);

            model.optimize();

            // for (int j = 0; j < instance.ncities; j++) {
            //     for (int k = 0; k < instance.ncities; k++) {
            //         if(vship[j][k].get(GRB.DoubleAttr.X) > 1e-6)
            //             System.out.println(j+"_"+k +vship[j][k].get(GRB.DoubleAttr.X));
            //     }
            // }

            // 一次迭代完成，需要更新次梯度 也即是 slack
            double [] slack = new double[instance.ncities];
            for (int j = 0; j < instance.ncities; j++) {
                for (int k = 0; k < instance.ncities; k++) {
                    slack[j] += vship[k][j].get(GRB.DoubleAttr.X);
                }
                slack[j] -= instance.demand[j];
            }

            // 判断更新下界，如果下界没有更新，就记录没有更新的次数
            // 注意要将目标函数中常数项重新添加
            if(model.get(GRB.DoubleAttr.ObjVal) + ud > LB+1e-6){
                LB = model.get(GRB.DoubleAttr.ObjVal)+ ud ;
                same = 0;
            }else {
                same++;
            }

            if(same == sameLimit){
                theta /= 2;
                same = 0;
            }

            // 计算sum(slack^2) ,即范数
            squareSum = 0;
            for (int j = 0; j < instance.ncities; j++) {
                squareSum += slack[j] * slack[j];
            }
            // 更新步长
            step = theta * (UB - (model.get(GRB.DoubleAttr.ObjVal)+ud)) / squareSum;

            // 更新梯度
            for (int j = 0; j < instance.ncities; j++) {
                if (xlambda[j] > step * slack[j]) {
                    // 往次梯度的反方向进行更新，如果其小于0，则令为0，因为lambda是非负的。
                    xlambda[j] -= step * slack[j];
                }else {
                    xlambda[j] = 0;
                }
            }

            selected_facility_supply = 0;
            for (int j = 0; j < instance.ncities; j++) {
                selected_facility_supply += vbuild[j].get(GRB.DoubleAttr.X) * instance.supply[j];
            }

            if(selected_facility_supply > sumDemand){
                // 说明当前迭代的方案是满足原始的条件的，因此恢复松弛的约束用于更新下界
                isLRModel = false;
                // 重新添加约束 sum(xij) >= dj for j in range(ncites)
                get_crelax();

                // 将模型中的yj变量的上下界设置为当前的解，即固定住
                for (int j = 0; j < instance.ncities; j++) {
                    vbuild[j].set(GRB.DoubleAttr.LB , vbuild[j].get(GRB.DoubleAttr.X));
                    vbuild[j].set(GRB.DoubleAttr.UB , vbuild[j].get(GRB.DoubleAttr.X));
                }

                model.setObjective(obj_shipCost , GRB.MINIMIZE);
                model.optimize();

                // 更新上界，如果当前的解有提升，就更新UB
                UB = Math.min(UB , model.get(GRB.DoubleAttr.ObjVal));

                // 重新设置yj的上下界
                for (int j = 0; j < instance.ncities; j++) {
                    vbuild[j].set(GRB.DoubleAttr.LB , 0);
                    vbuild[j].set(GRB.DoubleAttr.UB , 1);
                }
            }

            UBList.add(UB);
            LBList.add(LB);
            stepList.add(step);
            thetaList.add(theta);
            gapList.add((UB-LB)/UB*100);
        }
    }

    private GRBLinExpr get_obj_shipCost() {
        GRBLinExpr obj_shipCost = new GRBLinExpr();
        for (int i = 0; i < instance.ncities; i++) {
            for (int j = 0; j < instance.ncities; j++) {
                obj_shipCost.addTerm(instance.shipCost[i][j],vship[i][j]);
            }
        }
        return obj_shipCost;
    }

    private void get_crelax() throws GRBException {
        for (int j = 0; j < instance.ncities; j++) {
            GRBLinExpr expr1 = new GRBLinExpr();
            for (int i = 0; i < instance.ncities; i++) {
                expr1.addTerm(1,vship[i][j]);
            }
            crleax[j] = model.addConstr(expr1 , GRB.GREATER_EQUAL , instance.demand[j] , "crleax" + j);
        }
    }
    public void report(){
        System.out.println("            *** summary report ***");
        System.out.println("   iter        LB         UB       step      theta     gap");
        for (int i = 0; i < UBList.size(); i++) {
            System.out.printf("%6d %10.2f %10.2f %10.6f %10.6f %10.5f  %n",i,LBList.get(i),UBList.get(i),stepList.get(i),thetaList.get(i) ,gapList.get(i)) ;
        }
    }
    private double relaxLB() throws GRBException {
        modelRelax = model.relax();

        modelRelax.optimize();

        return modelRelax.get(GRB.DoubleAttr.ObjVal);
    }
    private void init() throws GRBException {
        GRBEnv env = new GRBEnv();
        model = new GRBModel(env);
        vship = new GRBVar[instance.ncities][instance.ncities];
        vbuild = new GRBVar[instance.ncities];
        csupply = new GRBConstr[instance.ncities];
        crleax = new GRBConstr[instance.ncities];
        for (int i = 0; i < instance.ncities; i++) {
            vbuild[i] = model.addVar(0,1,0,GRB.BINARY,"vbuild" + i);
            for (int j = 0; j < instance.ncities; j++) {
                vship[i][j] = model.addVar(0,instance.demand[j],0,GRB.INTEGER,"vship" + i + "_" + j);
            }
        }
        // 添加约束，sum(xij) <= supply[j] * yj for i in range(ncites)
        for (int i = 0; i < instance.ncities; i++) {
            GRBLinExpr expr = new GRBLinExpr();
            for (int j = 0; j < instance.ncities; j++) {
                expr.addTerm(1,vship[i][j]);
            }
            expr.addTerm(-1*instance.supply[i],vbuild[i]);
            csupply[i] = model.addConstr(expr , GRB.LESS_EQUAL , 0 , "csupply" + i);
        }

        GRBLinExpr expr = new GRBLinExpr();
        for (int i = 0; i < instance.ncities; i++) {
            expr.addTerm(1,vbuild[i]);
        }
        cbuildlimit = model.addConstr(expr , GRB.LESS_EQUAL , instance.buildlimit , "cbuildlimit");


        get_crelax();

        model.setObjective(get_obj_shipCost() , GRB.MINIMIZE);

        model.set(GRB.IntParam.OutputFlag , 0);

        model.update();
    }


}
