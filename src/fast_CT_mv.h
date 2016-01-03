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

        double  Get_diff(int,int );	//��J�y�ШD�t��
        double  Get_MVdiff(int,int); //��J�V�q�D�t��
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
   ����B�J
   Analysis_data()  Ū���
   Cell_analysis_data_64()  �t��CU size
   Enhance_cell() 64*64��m�[�j
   */
class MVmap_CT{
    public:
        MVmap_CT ();
        // create / destroy / initialize / copy    

        // member functions for Ū�J�B�]��
        void Count_frame();

        void Analysis_data();
        void Cell_analysis_data_64(int,int);
        void Enhance_cell(int,int,int,int);

        //�ϥ�MVmap
        int  Get_currframe(); //return _currframe
        int  Get_ReducedSR(int,int,int);       
        int  Get_ReducedSR_64(int,int,int);
        int  Get_MVnum_64(int,int);
        MV_CT  Get_nearestMV_64(int,int,int);
        //

        void creat_dec264pred(int,int,int,TComMv &); //�S�b�ΤF
        int  count ;
        int  min_diff;
    private:
        /*
           �o��S���ʺA�t�m
           _MVmap[frame to be encode]
           _MVmap[total frame * CUs per frame]
           �o��error�i�H�]�j�@�I�]
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






