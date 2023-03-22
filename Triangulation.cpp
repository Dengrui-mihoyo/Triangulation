#include<iostream>
#include<algorithm>
#include<cmath>
#include<vector>
using namespace std;

typedef struct                      //三维点的横纵列坐标
{                 
    double x;                           
    double y;
    double z;
    int index;                      //索引号
}Point;

typedef struct                      //三维点之间的边
{
    /* data */
    Point left;
    Point right;
    int count;                      //用于控制边的计数，若为0，则删除此边
}Edge;

typedef struct{
    //圆心的坐标点
    Point V[3];                     //三角形的三个顶点
    Edge E[3];                      //三角形的三条边
    double x_c;                     //外接圆的圆心的x坐标
    double y_c;                     //外接圆的圆心的y坐标
    double r;                       //外接圆的圆心的半径


}Triangle;




/*
三角剖分类
*/
class DelaunayTriangulation{
    public:
    vector<Point> point;                            //生成point数组
    vector<Triangle> triangle;                      //生成triangle数组
    vector<Edge> edge;                              
    vector<int> badTriangle;                   //不满足空圆性质的三角形


    DelaunayTriangulation(Point p1, Point p2, Point p3, Point p4);
    //构造初始大三角
    DelaunayTriangulation(){}                       //无参构造函数
    ~DelaunayTriangulation();                       //析构函数
    void Circumscribe(double &x_center,             //外接圆计算函数
                      double &y_center, 
                      double &r, 
                      Point v1, 
                      Point v2, 
                      Point v3);
    // void AddPoint(vector<Point> point);
    bool AddPoint(double x, double y, double z);
    void GenerateTri(Point v1, Point v2, Point v3);
    void DeleteTri(int n, vector<Triangle> triangle);
    bool Iscirle(double x, double y, Triangle triangle);
    void DeleteFrame();
    void Boundary_re(int startPoint, int endPoint);
    // bool operator ==(const Point &p)
    // {
    //     return (this->x == p.x ) && (this->y == p.y) && (this->z == p.z);
    // }


};
/***
 * @brief 构造函数构建矩形包围盒，生成两个超三角形作为初始三角网络
*/
DelaunayTriangulation::DelaunayTriangulation(Point p1, Point p2, Point p3, Point p4){
    // Point p_min = points[0];
    // Point p_max = points[0];
    // for (const auto & p :points){                       //遍历points上的节点
    //     if (p.x < p_min.x)  p_min.x = p.x;
    //     if (p.y < p_min.y)  p_min.y = p.y;
    //     if (p.z < p_min.z)  p_min.z = p.z;
    //     if (p.x > p_max.x)  p_max.x = p.x;
    //     if (p.y > p_max.y)  p_max.y = p.y;
    //     if (p.z > p_max.z)  p_max.z = p.z;
    //     找出最边界的顶点设定为最小和最大的坐标点
    //     GenerateTri(p1, p2, p3, p4);
    point.resize(4);
    edge.resize(4);
    point[0] = p1;
    point[1] = p2;
    point[2] = p3;
    point[3] = p4;
    Edge e1 = {0, 1, -1};             //描述顶点之间的边，此时并不存在
    Edge e2 = {1, 2, -1};
    Edge e3 = {2, 3, -1};
    Edge e4 = {0, 1, -1};                 
    edge[0] = e1;                     //将边的信息放入Edge中
    edge[1] = e2;
    edge[2] = e3;
    edge[3] = e4;
    GenerateTri(p1, p2, p3);          //生成两个超三角形作为初始三角网络
    GenerateTri(p1, p3, p4);
}


/***
 * @brief 析构函数实现
 *        将生成的vector数组长度都置为0，释放空间
*/
DelaunayTriangulation::~DelaunayTriangulation(){
    point.resize(0);
    triangle.resize(0);
    edge.resize(0);
}


/***
 * @brief Circumscribe实现:三角形外接圆的圆心的横纵坐标计算和半径计算
*/
void DelaunayTriangulation::Circumscribe(double &x_center, double &y_center, double &radius,  Point v1, Point v2,  Point v3){
    double x_1, x_2, x_3, y_1, y_2,y_3;
    //从点v1，v2，v3中获得横纵坐标
    x_1 = v1.x;
    y_1 = v1.y;
    x_2 = v2.x;
    y_2 = v2.y;
    x_3 = v3.x;
    y_3 = v3.y;
    /*具体公式在论文中体现*/
    x_center = (((y_2 - y_1)*(y_3*y_3 - y_1*y_1 + x_3*x_3 - x_1*x_1))
               -((y_3 - y_1)*(y_2*y_2 - y_1*y_1 + x_2*x_2 - x_1*x_1)))
               /((2*(x_3 - x_1)*(y_2 - y_1)) - (2*(x_2 - x_1)*(y_3 - y_1)));

    y_center = (((x_2 - x_1)*(x_3*x_3 - x_1*x_1 + y_3*y_3 - y_1*y_1))
               -((x_3 - x_1)*(x_2*x_2 - x_1*x_1 + y_2*y_2 - y_1*y_1)))
               /((2*(y_3 - y_1)*(x_2 - x_1)) - (2*(y_2 - y_1)*(x_3 - x_1)));

    radius = sqrt(pow(x_center - x_1, 2) + pow(y_center - y_1, 2));

}

/**
 * @brief GenerateTri函数:生成新的三角形函数实现
 *                        在三个顶点的外接圆内构造三角形
 */
void DelaunayTriangulation::GenerateTri(Point v1, Point v2, Point v3){
    double x_center;
    double y_center;
    double radius;
    Circumscribe(x_center, y_center, radius, v1, v2, v3);
    Triangle newtriangle = {(v1, v2, v3),
                            ((v1, v2, 1), (v1, v3, 1), (v2, v3, 1)), 
                            x_center, y_center, radius};
    triangle.push_back(newtriangle);
    int n = (int)edge.size();
    bool flag = true;                                                                //设置标记位（代表边是否在边集合中）
    for (int i = 0; i <= 2; i ++){
        for (int j = 0; j < n; j++){
            if (newtriangle.E[i].left.x == edge[j].left.x &&                         //也可以重载运算符
                newtriangle.E[i].left.y == edge[j].left.y &&
                newtriangle.E[i].left.z == edge[j].left.z &&
                newtriangle.E[i].right.x == edge[j].right.x &&
                newtriangle.E[i].right.x == edge[j].right.x &&
                newtriangle.E[i].right.x == edge[j].right.x &&
                edge[j].count != -1){
                    flag = false;
                    edge[j].count++;
                    break;
                }
            else if(newtriangle.E[i].left.x == edge[j].left.x &&
                    newtriangle.E[i].left.y == edge[j].left.y &&
                    newtriangle.E[i].left.z == edge[j].left.z &&
                    newtriangle.E[i].right.x == edge[j].right.x &&
                    newtriangle.E[i].right.x == edge[j].right.x &&
                    newtriangle.E[i].right.x == edge[j].right.x &&
                    edge[j].count == -1){
                        flag = false;
                        break;
                }
        }
        if (flag == 1){                                                              //如果flag没有变化的话，即新三角形的边不存在于Edge中
            edge.push_back(newtriangle.E[i]);                                        //需要将新三角的边加入Edege的Vector数组中
        }
    }
}

/***
 * @brief: AddPoint方法实现，用于加入点
 *         并且删除不合符条件的三角形，保留边框生成新的三角形
*/
bool DelaunayTriangulation::AddPoint(double x, double y, double z){
    vector<Edge> edgeFrame;                                                          //生成（删除三角形留下的）边框容器
    Point newPoint = {x,y,z};
    int n = (int)triangle.size();                                                    //badTriangle的索引号代表了坏三角的索引，若索引号< 0
    point.push_back(newPoint);                                                       //则代表坏三角形已经被删除point.push_back(newPoint);
    for (int i = 0; i < n ; i++){
        if (Iscirle(x,y,triangle[i]))
            badTriangle.push_back(i);
    }
    int mutex_1 = (int)badTriangle.size();
    for (int i = 0; i < mutex_1; i++){                                               //删除坏三角，但保留其边框，用于生成新的三角形
        // DeleteTri()
        /***
         * 未完成功能
        */
        for (int j = i + 1; j < mutex_1; j++){                                       //递增删除编号，保留0编号的三角形，生成新的三角形
            badTriangle[j] -= 1;
        }
    }
    int pointSize = (int)point.size();
    for (int i = 0; i < (int)edgeFrame.size(); i++){
        if (pointSize - 1 < edgeFrame[i].left.index){
            GenerateTri(point[pointSize - 1], 
                        point[edgeFrame[i].left.index], 
                        point[edgeFrame[i].right.index]);
        }
        //检查最大节点编号的和边框左右节点的大小，如果right > max, left > max则使用GenerateTri(left,right,max),left和right的位置可以互换
        else if(pointSize - 1 > edgeFrame[i].left.index && pointSize - 1 < edgeFrame[i].right.index){
            GenerateTri(point[edgeFrame[i].left.index], 
                        point[pointSize - 1], 
                        point[edgeFrame[i].right.index]);
        }
        else{
            GenerateTri(point[edgeFrame[i].left.index], 
                        point[edgeFrame[i].right.index], 
                        point[pointSize - 1]);
        }

    }
    if ((int)triangle.size() > n){
        return true;                    //添加完成
        cout<<"三角形添加成功，并加入了新的节点"<<endl;
    }
    else{
        return false;                   //添加失败
        cout<<"新的三角形生成失败"<<endl;
    }
}


/***
 * @brief: DeleteTri方法实现，用于删除坏三角形(不符合空圆性质的三角形)
*/
void DelaunayTriangulation::DeleteTri(int n, vector<Triangle> triangle){

}



/***
 * @brief: Iscircle方法实现，用于判断三角形的外接圆是否满足空圆性质
 *         若返回false则不在圆内，即满足空圆性质
 *         若返回True则在园内，不满足空圆性质
*/
bool DelaunayTriangulation::Iscirle(double x, double y, Triangle triangle){
    int distance = sqrt(pow(triangle.x_c - x,2) + pow(triangle.y_c - y,2));
    return distance < triangle.r;
}



/***
 * @brief: DeleteFrame方法实现，用于删除初始生成的四边形
*/
void DelaunayTriangulation::DeleteFrame(){
    vector<Edge> edgeFrame;
    //删除初始的四个顶点
    for (int i = 0; i < 4; i++){
        point.erase(point.begin());
    }
    //检查每一个三角形的第一个节点的索引值是否位0,1,2,3
    //如果是则删除其三角形并删除边框，让i的值--(三角形的个数减少)
    for (int i = 0; i < (int)triangle.size(); i++){
        if (triangle[i].V[0].index == 0 || 
            triangle[i].V[0].index == 1 || 
            triangle[i].V[0].index == 2 || 
            triangle[i].V[0].index == 3)
        {
            //DeleteTri
            /**
             * @brief 删除三角形功能还未规划实现
             */
            edgeFrame.resize(0);
            i--;
        }
        //将其他点的位置向前移4个
        else{
            for (int j = 0; j< 3; j++)
            {
                triangle[i].V[j].index -= 4;
                triangle[i].E[j].left.index -= 4;
                triangle[i].E[j].right.index -= 4;
            }
        }
    }
    for (int i =0; i < 4; i++){                                 //利用迭代器删除前四条边元素
        edge.erase(edge.begin());
    }
    int n = (int)edge.size();
    for (int i = 0; i < n; i++){                                //循环将每条边左右节点的索引值进行更改（删除了前四个节点）
        edge[i].left.index -= 4;
        edge[i].right.index -= 4;
    }
}


/**
 * @brief boundary_re()恢复边界函数实现
 *        用于实现超三角形的边界恢复
 */

void DelaunayTriangulation::Boundary_re(int startPoint, int endPoint){
    vector<Edge> edgeFrame;
    int n = edge.size();
    for (int i = 0; i < n; i++){
        if (triangle[i].V[0].index >= (startPoint-1) && triangle[i].V[2].index <= endPoint){
            /**
             * @brief deleteTri函数尚未实现
             * 
             */
            // DeleteTri()
            edgeFrame.resize(0);
            i--;
        }
    }
}