#include <string>
#include <vector>
#include <functional>
using namespace std;

//node 存MB input資訊
class Node
{
    public:
        Node();
        Node(int partition,int flag_intra,int flag_skip);
        //set
        void set_weight(int weight) {uweight=weight;};
        void set_depth(int depth)	{udepth=depth;};
        void set_part(int partition)	{partition_num=partition;};
        void set_intra(int intra)	{is_intra=intra;};
        void set_skip(int skip)	{is_skip=skip;};
        void set_equalMV(int MV){is_equalMV=MV;};
        void set_decMV(int MVx,int MVy){dec_MV[0]=MVx;dec_MV[1]=MVy;};
        void set_coff(double cof)	{coff=cof;};
        //get
        int get_partnum()	{return partition_num;};
        int get_weight()	{return uweight;};
        int get_depth()		{return udepth;};
        int isintra()		{return is_intra;};
        int isskip()		{return is_skip;};
        int isequalMV()      {return is_equalMV;};
        double get_coff()		{return coff;};
        int  MBadd;
        bool MB_activity;
        int  dec_MV [2]; //mode 用 dec MV ,0 = x,1 = y,
    private:
        double coff;
        int  is_intra;		//flag of intra
        int  is_skip;		//flag of skip
        int  is_equalMV;

        int  partition_num; //partition number of input MB
        int  uweight;		//total weight of this and under node (from partition)
        int  udepth;

};
//為了用find_if
class   Node_MBadd_equ   :public   unary_function <Node,bool> 
{ 
    public: 
        int   i; 
        explicit   Node_MBadd_equ(const   int   in):   i(in){} 
        bool     operator   ()(const   Node&   node){return   node.MBadd==i;} 
}; 


//pre CU scale
class MODEtree_CT
{
    public:
        MODEtree_CT();
        //初始、刪除
        void create();
        void destroy();

        //讀入、處理
        void set_node(int uiIdx,Node inNode){dec264_tree[uiIdx]=inNode;};
        void set_node_false();
        //	void Analysis_data();

        //使用
        Node get_node(int uiIdx){return dec264_tree[uiIdx];};
        void show_node(int uiIdx);
        int    count_SubMB;
        int    equalMV;
        int    equalMV_d1[4];
    private:
        Node*  dec264_tree; // array of dec264 tree




};

//pre frame scale
class MODEframe_tree_CT
{
    public:
        MODEframe_tree_CT();
        //初始、刪除
        void create(int CU_perFrame);
        void destroy();

        //讀入、處理
        void set_tree_dec264(std::vector<Node> dec264_node,int CUadd,int MBmap[16]);
        void set_tree_dec264(int CUadd);
        MODEtree_CT get_tree(int uiIdx){return dec264_CU_tree[uiIdx];};
        //使用

    private:
        void count_frameMBActivity();

        MODEtree_CT*  dec264_CU_tree; // array of dec264 tree
        int  _CU_perFrame;
};






class Transcoder_CT_mode2
{
    public:
        Transcoder_CT_mode2();
        void init(int width,int height,int FrameToEncode,int QP,int flag_ZB);
        //初始、刪除
        void create();
        void destroy();

        //處理
        void  set_Tcoder_CU_AMBA(int currframe,int currCUadd);
        int  dec264_isskip[16];   
        int  dec264_part[16];
        int  dec264_subMB;
        int  dec264_skip;
        int  dec264_totalWeight;
        double  dec264_AMBA_d0;  //Avg MBA using in depth 0
        double  dec264_AMBA_d1[4];//Avg MBA using in depth 1
        double  dec264_AMBA_d2[16];//Avg MBA using in depth 2

        int  dec264_ZEROMV_d0;//0 MV using in depth 1
        int  dec264_ZEROMV_d1[4];//0 MV using in depth 2

        //使用
        double get_frameMBA(int frame){return avg_MBA[frame];}
        MODEtree_CT get_tree(int currframe,int currCUadd);
        void show(int frame);
        void MBtoCU_mapping(const int CUadd,int (&MBmap)[16]); //
        int get_QP(){return _QP;};
        int _width,_height,_totalFrame;//基本資訊
    private:
        void Analysis_data();//從檔案讀入資料
        void Processing_data_CU(std::vector <Node> dec264_node);//建立CU tree
        MODEframe_tree_CT* _currFrame;//存某frame裡面的資訊
        int _Count_currframe;
        int _flag_ZB,_QP;
        int RowCell_MB,ColCell_MB,FrameShift_MB;//算出一個frame有幾個MB
        int RowCell_CU,ColCell_CU,FrameShift_CU;//算出一個frame有幾個CU
        double total_MBA;
        double* avg_MBA;		

};
