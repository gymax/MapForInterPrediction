#include "../../Lib/TlibCommon/CommonDef.h"
#include <string>
#include <vector>
#include "../../Lib/TLibCommon/TComMv.h"
#include "../../Lib/TlibCommon/CommonDef.h"
using namespace std;
class MV_CT {
    public:
        MV_CT();
        MV_CT(int,int,int,int);
        int		add_X,add_Y;
        int		vector_X,vector_Y;

        double  Get_diff(int,int );	//輸入座標求差值
        double  Get_MVdiff(int,int); //輸入向量求差值
        double  Get_radius(); 

    private:

        double _diff;
        double _radius;
        /*
           ---------------operator-----------------------
           */
    public:
        const MV_CT& operator += (const MV_CT& rcMv)
        {
            vector_X += rcMv.vector_X;
            vector_X += rcMv.vector_X;
            return  *this;
        }

        const MV_CT& operator-= (const MV_CT& rcMv)
        {
            vector_X -= rcMv.vector_X;
            vector_Y -= rcMv.vector_Y;
            return  *this;
        }

        const MV_CT& operator>>= (const Int i)
        {
            vector_X >>= i;
            vector_Y >>= i;
            return  *this;
        }

        const MV_CT& operator<<= (const Int i)
        {
            vector_X <<= i;
            vector_Y <<= i;
            return  *this;
        }


        Bool operator== ( const MV_CT& rcMv ) const
        {
            return (vector_X==rcMv.vector_X && vector_Y==rcMv.vector_Y);
        }

        Bool operator!= ( const MV_CT& rcMv ) const
        {
            return (vector_X!=rcMv.vector_X || vector_Y!=rcMv.vector_Y);
        }
}; 

/*
   執行步驟
   Analysis_data()  讀資料
   Cell_analysis_data_64()  配成CU size
   Enhance_cell() 64*64位置加強
   */
class MVmap_CT{
    public:
        MVmap_CT ();
        // create / destroy / initialize / copy    

        // member functions for 讀入、設值
        void Count_frame();

        void Analysis_data();
        void Cell_analysis_data_64(int,int);
        void Enhance_cell(int,int,int,int);

        //使用MVmap
        int  Get_currframe(); //return _currframe
        int  Get_ReducedSR(int,int,int);       
        int  Get_ReducedSR_64(int,int,int);
        int  Get_MVnum_64(int,int);
        MV_CT  Get_nearestMV_64(int,int,int);
        //

        void creat_dec264pred(int,int,int,TComMv &); //沒在用了
        int  count ;
        int  min_diff;
    private:
        /*
           這邊沒有動態配置
           _MVmap[frame to be encode]
           _MVmap[total frame * CUs per frame]
           發生error可以設大一點跑
           */
        MV_CT		   _currMV;
        vector<MV_CT>  _MVmap[110];
        vector<MV_CT>  _MVmap_64[90000];
        int	 _currframe;	
        int  _inputframe;
        //int  _frameToEnco;
        int RowCell,ColCell,FrameShift;
        void CU64mv_average(vector<MV_CT> &CUmv,int All_CUadd);
        void CU64mv_median(vector<MV_CT> &CUmv,int All_CUadd);
        void CU64mv_mode(vector<MV_CT> &CUmv,int All_CUadd);
        void CU32mv_average(vector<MV_CT> &CUmv,int All_CUadd);
        void CU32mv_median(vector<MV_CT> &CUmv,int All_CUadd);
        void CU32mv_mode(vector<MV_CT> &CUmv,int All_CUadd);
};






