#include <iostream> 
#include <fstream>
#include <vector>
#include <algorithm>
#include "fast_CT_mv.h" 
#include <math.h>
//block offset
/*  order       ________
    | 1  2 |
    | 3  4 |
    --------		*/
//#define  offset32X32  16
//#define  offset64X64  32
#define  offset64       64
#define  offset64xy     32
#define  offsetMB       16
static	int offset32X32x[] = {16,48,16,48};
static  int offset32X32y[] = {16,16,48,48};


//MV的架構成員
MV_CT::MV_CT(){
    _radius  =  0.0;
    add_X=0;
    add_Y=0;
    vector_X=0;
    vector_Y=0;
}

MV_CT::MV_CT(int addx, int addy, int vx, int vy){
    add_X=addx;
    add_Y=addy;
    vector_X=vx;
    vector_Y=vy;
    _radius=Get_radius(); 
}

double MV_CT::Get_diff(int input_X,int input_Y){

    _diff=ceil(sqrt(pow(double(add_X-input_X),2)+pow(double(add_Y-input_Y),2)));
    return _diff; 
}

double MV_CT::Get_MVdiff(int input_X,int input_Y){
    double MVdiff;
    MVdiff=ceil(sqrt(pow(double(vector_X-input_X),2)+pow(double(vector_Y-input_Y),2)));
    return MVdiff; 
}

double MV_CT::Get_radius(){

    _radius=ceil(sqrt(pow(double(vector_X),2)+pow(double(vector_Y),2)));
    return _radius;
}

//MVMAP 函式
MVmap_CT::MVmap_CT(){
    _currframe	=	-1;
    count		=	0;
}



void MVmap_CT::Count_frame(){
    _currframe++;
}

int MVmap_CT::Get_currframe(){
    return _currframe;
}
void MVmap_CT::Analysis_data(){
    ifstream fin("dec_264_mv.txt");
    string str;
    int currframe	=	0;
    int selectframe	=	0;
    int frame_end	=   0;
    //對frame 生成_MVmap
    //讀檔

    for(int i = 0 ;i<101 ; i++)
    {
        while(fin.good())
        {

            fin >> str;
            //判斷frame是否結束
            if((str[0]=='-')&(frame_end>0))
            {	
                break;  		  
            }

            //判斷是第幾個frame
            if(str[0]=='{')
            { 
                fin >> str;
                frame_end = 0;
                currframe++;
                frame_end++;
                //cout << "目前是frame : " << currframe << endl;
            }

            //read  MV & ADD


            if(str[0]=='<')
            {
                string  strMVx,strMVy;

                fin >> strMVx;
                fin >> strMVy;

                _currMV.vector_X=atoi(strMVx.c_str());
                _currMV.vector_Y=atoi(strMVy.c_str());
                //MVx=atoi(strMVx.c_str());
                //MVy=atoi(strMVy.c_str());

                //cout <<"mv("<<_currMV.vector_X<<","<<_currMV.vector_Y<<")";

            }

            if(str[0]=='(')
            {
                string  strADDx,strADDy;

                fin >> strADDx;
                fin >> strADDy;

                _currMV.add_X=atoi(strADDx.c_str());
                _currMV.add_Y=atoi(strADDy.c_str());
                _MVmap[currframe].push_back(_currMV);
            }

        }

        _inputframe=currframe;
    }
    fin.close();
};


int MVmap_CT::Get_ReducedSR(int selectframe,int input_X,int input_Y){

    //計算最近點及MV radius
    //將距離存入diffADD 
    int		diff;
    int   MIN_VAULE	= 8000;
    int radius,cloest_X,cloest_Y;

    for(vector<MV_CT>::iterator it = _MVmap[selectframe].begin(); it != _MVmap[selectframe].end(); it++) 
    {

        MV_CT	getMV=*it;
        diff=ceil(sqrt(pow(double(getMV.add_X-input_X),2)+pow(double(getMV.add_Y-input_Y),2)));

        //找出最小值
        if (diff<MIN_VAULE)
        {	
            radius=ceil(sqrt(pow(double(getMV.vector_X),2)+pow(double(getMV.vector_Y),2)));
            cloest_X=getMV.add_X;
            cloest_Y=getMV.add_Y;
            MIN_VAULE=diff;
        }	
        if(MIN_VAULE <= 1)break;
    }     
    //cout<<"離("<<input_X<<","<<input_Y<<")最近的點是"<<"("<<cloest_X<<","<<cloest_Y<<")"<<endl;
    min_diff=MIN_VAULE;

    return radius;
}


void MVmap_CT::creat_dec264pred(int selectframe,int input_X,int input_Y,TComMv&  rcMvPred){

    //回傳264 dec MVP 
    int		diff;
    int   MIN_VAULE	= 8000;
    int   cloest_X,cloest_Y;

    for(vector<MV_CT>::iterator it = _MVmap[selectframe].begin(); it != _MVmap[selectframe].end(); it++) 
    {

        MV_CT	getMV=*it;
        diff=ceil(sqrt(pow(double(getMV.add_X-input_X),2)+pow(double(getMV.add_Y-input_Y),2)));

        //找出最小值
        if (diff<MIN_VAULE)
        {	
            cloest_X=getMV.vector_X;
            cloest_Y=getMV.vector_X;
            MIN_VAULE=diff;
        }	
        if(MIN_VAULE <= 1)break;
    }
    rcMvPred.set(short(cloest_X)/4,short(cloest_Y)/4);
}

void MVmap_CT::Cell_analysis_data_64(int width,int height){

    RowCell=ceil(double(width/64.0));
    ColCell=ceil(double(height/64.0));
    FrameShift=RowCell*ColCell;
    for(int frame=0 ; frame < _inputframe ; frame++){
        // frame+1 是一個shift 因為MVmap[0]是i frame
        for(vector<MV_CT>::iterator it = _MVmap[frame+1].begin(); it != _MVmap[frame+1].end(); it++) 
        {  
            MV_CT	getMV=*it;
            int RowStep,ColStep,CUadd,All_CUadd;
            RowStep=int(getMV.add_X/64);
            ColStep=int(getMV.add_Y/64);
            CUadd=ColStep*RowCell+RowStep;
            All_CUadd=CUadd+frame*FrameShift;
            _MVmap_64[All_CUadd].push_back(getMV);
        }
    }
}

int MVmap_CT::Get_ReducedSR_64(int selectframe,int input_X,int input_Y){

    //計算最近點及MV radius
    //將距離存入diffADD 
    int		diff;
    int   MIN_VAULE	= 8000;
    int radius;
    int RowStep,ColStep,CUadd,All_CUadd;

    RowStep=int(input_X/64);
    ColStep=int(input_Y/64);
    CUadd=ColStep*RowCell+RowStep;
    //selectframe-1是因為要對齊
    All_CUadd=CUadd+(selectframe-1)*FrameShift;


    for(vector<MV_CT>::iterator it = _MVmap_64[All_CUadd].begin(); it != _MVmap_64[All_CUadd].end(); it++) 
    {
        MV_CT   bestMV;
        MV_CT	getMV=*it;
        diff=getMV.Get_diff(input_X,input_Y);

        //找出最小值
        if (diff<MIN_VAULE)
        {	
            radius=getMV.Get_radius();
            MIN_VAULE=diff;
        }	
        if(MIN_VAULE <= 1)break;
    }     

    min_diff=MIN_VAULE;

    return radius;
}

MV_CT MVmap_CT::Get_nearestMV_64(int selectframe,int input_X,int input_Y){

    //計算最近點及MV radius
    //將距離存入diffADD 
    int		diff;
    int   MIN_VAULE	= 8000;
    int radius;
    int RowStep,ColStep,CUadd,All_CUadd;
    MV_CT   bestMV;

    RowStep=int(input_X/64);
    ColStep=int(input_Y/64);
    CUadd=ColStep*RowCell+RowStep;
    //selectframe-1是因為要對齊
    All_CUadd=CUadd+(selectframe-1)*FrameShift;

    for(vector<MV_CT>::iterator it = _MVmap_64[All_CUadd].begin(); it != _MVmap_64[All_CUadd].end(); it++) 
    {

        MV_CT	getMV=*it;

        diff=getMV.Get_diff(input_X,input_Y);

        //找出最小值
        if (diff<MIN_VAULE)
        {	
            radius=getMV.Get_radius();
            MIN_VAULE=diff;
            bestMV=getMV;
        }	
        if(MIN_VAULE <= 1)break;
    }     

    min_diff=MIN_VAULE;

    return bestMV;
}


int MVmap_CT::Get_MVnum_64(int selectframe,int CUadd)
{

    int All_CUadd=	CUadd+(selectframe-1)*FrameShift; 

    if(selectframe>0){
        return _MVmap_64[All_CUadd].size();
    }else{
        return 0;
    }
}
/*
Enhance_cell
提供32x32 64x64更好的預測
method 0  -> average
method 1  -> 中位數
method 2  -> 眾數
*/

void MVmap_CT::Enhance_cell(int width, int height,int method_64,int method_32)
{  


    for(int frame=0 ; frame < (_inputframe-1) ; frame++)
        //_inputframe-1是因為要扣掉i frame
    {
        for(int CUadd =0; CUadd <FrameShift; CUadd++)
        {

            int All_CUadd=CUadd+frame*FrameShift;
            if(!_MVmap_64[All_CUadd].empty())
            {
                switch(method_64)
                {//對CU 64x64 做處理
                    case 0:
                        CU64mv_average(_MVmap_64[All_CUadd],All_CUadd);  

                        break;

                    case 1:
                        CU64mv_median(_MVmap_64[All_CUadd],All_CUadd);  //median 有做Nx2N 2NxN 點
                        break;

                    case 2:
                        //CU64mv_mode(_MVmap_64[All_CUadd],All_CUadd);
                        break;
                }

                switch(method_32)
                {//對CU 32x32 做處理
                    case 0:
                        CU32mv_average(_MVmap_64[All_CUadd],All_CUadd);  
                        break;

                    case 1:
                        CU32mv_median(_MVmap_64[All_CUadd],All_CUadd);  //median 有做Nx2N 2NxN 點
                        break;

                    case 2:
                        //CU32mv_mode(_MVmap_64[All_CUadd],All_CUadd);
                        break;
                }
            }
        }//CU

    }//frame
}


void MVmap_CT::CU64mv_average (vector<MV_CT> &CUmv,int All_CUadd){
    int CU_offset_x = 0;
    int CU_offset_y = 0;
    int currCU_addx;
    int currCU_addy;
    int currCUadd=All_CUadd %FrameShift;
    vector<int> MBmvx;
    vector<int> MBmvy;
    CU_offset_x=(currCUadd)%(RowCell);
    CU_offset_y=(currCUadd-(currCUadd)%(RowCell))/RowCell;
    currCU_addx=CU_offset_x*offset64;
    currCU_addy=CU_offset_y*offset64;


    /* frame size
       for(vector<MV_CT>::iterator it = CUmv.begin(); it != CUmv.end(); it++) 
       {
       MV_CT	getMV=*it;
       MBmvx.push_back(getMV.vector_X);
       MBmvy.push_back(getMV.vector_Y);
       }
       */

    for(int MBcount=0;MBcount<15;MBcount++)
    {
        int MB_offset_x=(MBcount)%(4);
        int MB_offset_y=(MBcount-(MBcount)%(4))/4;
        int     total_MBMVx=0,total_MBMVy=0;
        int     MVcount=0;
        for(vector<MV_CT>::iterator it = CUmv.begin(); it != CUmv.end(); it++) 
        {
            MV_CT	getMV=*it;

            if((getMV.add_X-currCU_addx)>(MB_offset_x*offsetMB)&&(getMV.add_X-currCU_addx)<((MB_offset_x+1)*offsetMB))
            {
                if((getMV.add_Y-currCU_addy)>(MB_offset_y*offsetMB)&&(getMV.add_Y-currCU_addy)<((MB_offset_y+1)*offsetMB))
                {
                    total_MBMVx+=getMV.vector_X;
                    total_MBMVy+=getMV.vector_Y;
                    MVcount++;
                }
            }
        }
        if(MVcount!=0)
        {
            MBmvx.push_back(total_MBMVx/MVcount);
            MBmvy.push_back(total_MBMVy/MVcount);
        }
    }

    int total_x=0,total_y=0;

    for(vector<int>::iterator it = MBmvx.begin(); it != MBmvx.end(); it++)
    {
        total_x+=*it;
    }
    for(vector<int>::iterator it = MBmvy.begin(); it != MBmvy.end(); it++)
    {
        total_y+=*it;
    }


    int n = MBmvx.size();
    int avg_x;
    int avg_y;
    if(n!=0)
    {
        avg_x=total_x/n;
        avg_y=total_y/n;
    }
    else
    {
        avg_x=0;
        avg_y=0;
    }

    MV_CT  enhanceMV_64x64(currCU_addx+offset64xy,currCU_addy+offset64xy,avg_x,avg_y);

    CUmv.push_back(enhanceMV_64x64);
}





void MVmap_CT::CU64mv_median (vector<MV_CT> &CUmv,int All_CUadd){
    int CU_offset_x = 0;
    int CU_offset_y = 0;
    int currCU_addx;
    int currCU_addy;
    int currCUadd=All_CUadd %FrameShift;
    vector<int> MBmvx;
    vector<int> MBmvy;
    CU_offset_x=(currCUadd)%(RowCell);
    CU_offset_y=(currCUadd-(currCUadd)%(RowCell))/RowCell;
    currCU_addx=CU_offset_x*offset64;
    currCU_addy=CU_offset_y*offset64;


    //先做 2Nx2N

    for(vector<MV_CT>::iterator it = CUmv.begin(); it != CUmv.end(); it++) 
    {
        MV_CT	getMV=*it;
        MBmvx.push_back(getMV.vector_X);
        MBmvy.push_back(getMV.vector_Y);
    }


    int n = MBmvx.size() / 2; 

    MV_CT  enhanceMV_64x64(currCU_addx+offset64xy,currCU_addy+offset64xy, 0, 0);
    nth_element(MBmvx.begin(), MBmvx.begin()+n, MBmvx.end()); 
    nth_element(MBmvy.begin(), MBmvy.begin()+n, MBmvy.end()); 
    enhanceMV_64x64.vector_X=MBmvx[n];
    enhanceMV_64x64.vector_Y=MBmvy[n];

    if(n!=0)  CUmv.push_back(enhanceMV_64x64);
#if fast_CT_mv_2NxN
    // 2NxN *2
    //清空 and 再宣告一組vector
    MBmvx.clear();
    MBmvy.clear();
    vector<int> MBmvx_2;
    vector<int> MBmvy_2;


    for(vector<MV_CT>::iterator it = CUmv.begin(); it != CUmv.end(); it++) 
    {
        MV_CT	getMV=*it;
        if(getMV.add_Y<currCU_addy+offset64xy)//2XxN:比y軸
        {//上面那群
            MBmvx.push_back(getMV.vector_X);
            MBmvy.push_back(getMV.vector_Y);
        }else{//下面那群
            MBmvx_2.push_back(getMV.vector_X);
            MBmvy_2.push_back(getMV.vector_Y);
        }
    }


    n = MBmvx. size() / 2; 
    int n_2 = MBmvx_2.size() / 2; //多令一個 n 
    if(n!=0)
    {
        MV_CT  enhanceMV_64x32_1(currCU_addx+32,currCU_addy+16, 0, 0);
        nth_element(MBmvx.begin(), MBmvx.begin()+n, MBmvx.end()); 
        nth_element(MBmvy.begin(), MBmvy.begin()+n, MBmvy.end()); 
        enhanceMV_64x32_1.vector_X=MBmvx[n];
        enhanceMV_64x32_1.vector_Y=MBmvy[n];
        CUmv.push_back(enhanceMV_64x32_1);
    }
    if(n_2!=0)
    {  
        MV_CT  enhanceMV_64x32_2(currCU_addx+32,currCU_addy+48, 0, 0);
        nth_element(MBmvx_2.begin(), MBmvx_2.begin()+n_2, MBmvx_2.end()); 
        nth_element(MBmvy_2.begin(), MBmvy_2.begin()+n_2, MBmvy_2.end());   
        enhanceMV_64x32_2.vector_X=MBmvx_2[n_2];
        enhanceMV_64x32_2.vector_Y=MBmvy_2[n_2]; 
        CUmv.push_back(enhanceMV_64x32_2);
    }

    // Nx2N *2
    //清空
    MBmvx.clear();
    MBmvy.clear();
    MBmvx_2.clear();
    MBmvy_2.clear();

    for(vector<MV_CT>::iterator it = CUmv.begin(); it != CUmv.end(); it++) 
    {
        MV_CT	getMV=*it;
        if(getMV.add_X<currCU_addx+offset64xy)//2XxN:比y軸
        {//左面那群
            MBmvx.push_back(getMV.vector_X);
            MBmvy.push_back(getMV.vector_Y);
        }else{//又面那群
            MBmvx_2.push_back(getMV.vector_X);
            MBmvy_2.push_back(getMV.vector_Y);
        }
    }


    n = MBmvx.size() / 2; 
    n_2 = MBmvx_2.size() / 2; //多令一個 n 
    if(n!=0)
    {


        MV_CT  enhanceMV_32x64_1(currCU_addx+16,currCU_addy+32, 0, 0);  
        nth_element(MBmvx.begin(), MBmvx.begin()+n, MBmvx.end()); 
        nth_element(MBmvy.begin(), MBmvy.begin()+n, MBmvy.end());   

        enhanceMV_32x64_1.vector_X=MBmvx[n];
        enhanceMV_32x64_1.vector_Y=MBmvy[n];  

        CUmv.push_back(enhanceMV_32x64_1);  
    }

    if(n_2!=0) {
        MV_CT  enhanceMV_32x64_2(currCU_addx+48,currCU_addy+32, 0, 0);
        nth_element(MBmvx_2.begin(), MBmvx_2.begin()+n_2, MBmvx_2.end()); 
        nth_element(MBmvy_2.begin(), MBmvy_2.begin()+n_2, MBmvy_2.end()); 


        enhanceMV_32x64_2.vector_X=MBmvx_2[n_2];
        enhanceMV_32x64_2.vector_Y=MBmvy_2[n_2];


        CUmv.push_back(enhanceMV_32x64_2);

    }
#endif


}


void MVmap_CT::CU64mv_mode (vector<MV_CT> &CUmv,int All_CUadd){
    int CU_offset_x = 0;
    int CU_offset_y = 0;
    int currCU_addx;
    int currCU_addy;
    int currCUadd=All_CUadd %FrameShift;
    vector<int> MBmvx;
    vector<int> MBmvy;
    CU_offset_x=(currCUadd)%(RowCell);
    CU_offset_y=(currCUadd-(currCUadd)%(RowCell))/RowCell;
    currCU_addx=CU_offset_x*offset64;
    currCU_addy=CU_offset_y*offset64;

    for(int MBcount=0;MBcount<15;MBcount++){
        int MB_offset_x=(MBcount)%(4);
        int MB_offset_y=(MBcount-(MBcount)%(4))/4;
        int     total_MBMVx=0,total_MBMVy=0;
        int     MVcount=0;
        for(vector<MV_CT>::iterator it = CUmv.begin(); it != CUmv.end(); it++) 
        {
            MV_CT	getMV=*it;

            if((getMV.add_X-currCU_addx)>(MB_offset_x*offsetMB)&&(getMV.add_X-currCU_addx)<((MB_offset_x+1)*offsetMB)){
                if((getMV.add_Y-currCU_addy)>(MB_offset_y*offsetMB)&&(getMV.add_Y-currCU_addy)<((MB_offset_y+1)*offsetMB)){
                    total_MBMVx+=getMV.vector_X;
                    total_MBMVy+=getMV.vector_Y;
                    MVcount++;
                }
            }
        }
        if(MVcount!=0){
            MBmvx.push_back(total_MBMVx/MVcount);
            MBmvy.push_back(total_MBMVy/MVcount);
        }
    }

    int n = MBmvx.size() / 2;    
    nth_element(MBmvx.begin(), MBmvx.begin()+n, MBmvx.end()); 
    nth_element(MBmvy.begin(), MBmvy.begin()+n, MBmvy.end()); 



    MV_CT  enhanceMV_64x64(currCU_addx+offset64xy,currCU_addy+offset64xy,0,0);

    CUmv.push_back(enhanceMV_64x64);
}




void MVmap_CT::CU32mv_average (vector<MV_CT> &CUmv,int All_CUadd){
    int CU_offset_x = 0;
    int CU_offset_y = 0;
    int currCU_addx;
    int currCU_addy;
    int currCUadd=All_CUadd %FrameShift;
    //用4個MBmv存CU 32*32
    vector<int> MBmvx[4];
    vector<int> MBmvy[4];
    //查表用的MB array mapping
    int MBmapping[16]={0,0,1,1,0,0,1,1,2,2,3,3,2,2,3,3};
    CU_offset_x=(currCUadd)%(RowCell);
    CU_offset_y=(currCUadd-(currCUadd)%(RowCell))/RowCell;
    currCU_addx=CU_offset_x*offset64;
    currCU_addy=CU_offset_y*offset64;


    /* frame size
       for(vector<MV_CT>::iterator it = CUmv.begin(); it != CUmv.end(); it++) 
       {
       MV_CT	getMV=*it;
       MBmvx.push_back(getMV.vector_X);
       MBmvy.push_back(getMV.vector_Y);
       }
       */

    for(int MBcount=0;MBcount<15;MBcount++){
        int MB_offset_x=(MBcount)%(4);
        int MB_offset_y=(MBcount-(MBcount)%(4))/4;
        int     total_MBMVx=0,total_MBMVy=0;
        int     MVcount=0;
        for(vector<MV_CT>::iterator it = CUmv.begin(); it != CUmv.end(); it++) 
        {
            MV_CT	getMV=*it;

            if((getMV.add_X-currCU_addx)>(MB_offset_x*offsetMB)&&(getMV.add_X-currCU_addx)<((MB_offset_x+1)*offsetMB)){
                if((getMV.add_Y-currCU_addy)>(MB_offset_y*offsetMB)&&(getMV.add_Y-currCU_addy)<((MB_offset_y+1)*offsetMB)){
                    total_MBMVx+=getMV.vector_X;
                    total_MBMVy+=getMV.vector_Y;
                    MVcount++;
                }
            }
        }	

        if(MVcount!=0){
            MBmvx[MBmapping[MBcount]].push_back(total_MBMVx/MVcount);
            MBmvy[MBmapping[MBcount]].push_back(total_MBMVy/MVcount);
        }
    }

    //CU中有4個32x32
    for(int CU32x32=0;CU32x32<4;CU32x32++){
        int total_x=0,total_y=0;

        for(vector<int>::iterator it = MBmvx[CU32x32].begin(); it != MBmvx[CU32x32].end(); it++){
            total_x+=*it;
        } 
        for(vector<int>::iterator it = MBmvy[CU32x32].begin(); it != MBmvy[CU32x32].end(); it++){
            total_y+=*it;
        } 
        int n = MBmvx[CU32x32].size();
        int avg_x;
        int avg_y;
        if(n!=0){
            avg_x=total_x/n;
            avg_y=total_y/n;
        }else{
            avg_x=0;
            avg_y=0;
        }

        MV_CT  enhanceMV_32x32(currCU_addx+offset32X32x[CU32x32],currCU_addy+offset32X32y[CU32x32],avg_x,avg_y);

        CUmv.push_back(enhanceMV_32x32);
    }
}


void MVmap_CT::CU32mv_median (vector<MV_CT> &CUmv,int All_CUadd){
    int CU_offset_x = 0;
    int CU_offset_y = 0;
    int currCU_addx;
    int currCU_addy;
    int currCUadd=All_CUadd %FrameShift;
    //用4個MBmv存CU 32*32
    vector<int> MBmvx[4];
    vector<int> MBmvy[4];
    //查表用的MB array mapping
    //int MBmapping[16]={0,0,1,1,0,0,1,1,2,2,3,3,2,2,3,3};
    CU_offset_x=(currCUadd)%(RowCell);
    CU_offset_y=(currCUadd-(currCUadd)%(RowCell))/RowCell;
    currCU_addx=CU_offset_x*offset64;
    currCU_addy=CU_offset_y*offset64;

    //32x32

    for(vector<MV_CT>::iterator it = CUmv.begin(); it != CUmv.end(); it++) 
    {
        MV_CT	getMV=*it;
        if((getMV.add_X-currCU_addx)<32&&(getMV.add_Y-currCU_addy)<32)
        { MBmvx[0].push_back(getMV.vector_X);
            MBmvy[0].push_back(getMV.vector_Y);}
        else if((getMV.add_X-currCU_addx)>32&&(getMV.add_Y-currCU_addy)>32)
        { MBmvx[3].push_back(getMV.vector_X);
            MBmvy[3].push_back(getMV.vector_Y);}
        else if((getMV.add_X-currCU_addx)>32)
        { MBmvx[1].push_back(getMV.vector_X);
            MBmvy[1].push_back(getMV.vector_Y);}
        else
        { MBmvx[2].push_back(getMV.vector_X);
            MBmvy[2].push_back(getMV.vector_Y);}
    }


    for(int CU32x32=0;CU32x32<4;CU32x32++)
    {
        int n = MBmvx[CU32x32].size() / 2;  
        MV_CT  enhanceMV_32x32(currCU_addx+offset32X32x[CU32x32],currCU_addy+offset32X32y[CU32x32], 0, 0);
        nth_element(MBmvx[CU32x32].begin(), MBmvx[CU32x32].begin()+n, MBmvx[CU32x32].end()); 
        nth_element(MBmvy[CU32x32].begin(), MBmvy[CU32x32].begin()+n, MBmvy[CU32x32].end()); 

        if(n!=0)
        {
            enhanceMV_32x32.vector_X=MBmvx[CU32x32][n];
            enhanceMV_32x32.vector_Y=MBmvy[CU32x32][n];
            CUmv.push_back(enhanceMV_32x32);
        }

    }//end for cu32
#if fast_CT_mv_2NxN
    //32x16
    //清空並再令一組vector
    for(int i=0;i<4;i++)
    {
        MBmvx[i].clear();
        MBmvy[i].clear();
    }
    vector<int> MBmvx_2[4];
    vector<int> MBmvy_2[4];


    for(vector<MV_CT>::iterator it = CUmv.begin(); it != CUmv.end(); it++) 
    {
        MV_CT	getMV=*it;
        if((getMV.add_X-currCU_addx)<32&&(getMV.add_Y-currCU_addy)<32)
        { 
            if(getMV.add_Y<currCU_addy+16)//比Y
            {//上面
                MBmvx[0].push_back(getMV.vector_X);
                MBmvy[0].push_back(getMV.vector_Y);
            }else{//下面
                MBmvx_2[0].push_back(getMV.vector_X);
                MBmvy_2[0].push_back(getMV.vector_Y);			
            }
        }
        else if((getMV.add_X-currCU_addx)>32&&(getMV.add_Y-currCU_addy)>32)
        { 
            if(getMV.add_Y<currCU_addy+48)
            {
                MBmvx[3].push_back(getMV.vector_X);
                MBmvy[3].push_back(getMV.vector_Y);
            }else{
                MBmvx_2[3].push_back(getMV.vector_X);
                MBmvy_2[3].push_back(getMV.vector_Y);			
            }
        }
        else if((getMV.add_X-currCU_addx)>32)
        { 
            if(getMV.add_Y<currCU_addy+16)
            {
                MBmvx[1].push_back(getMV.vector_X);
                MBmvy[1].push_back(getMV.vector_Y);
            }else{
                MBmvx_2[1].push_back(getMV.vector_X);
                MBmvy_2[1].push_back(getMV.vector_Y);			
            }
        }
        else
        { 
            if(getMV.add_Y<currCU_addy+48)
            {
                MBmvx[2].push_back(getMV.vector_X);
                MBmvy[2].push_back(getMV.vector_Y);
            }else{
                MBmvx_2[2].push_back(getMV.vector_X);
                MBmvy_2[2].push_back(getMV.vector_Y);			
            }
        }
    }


    for(int CU32x32=0;CU32x32<4;CU32x32++)
    {
        int n = MBmvx[CU32x32].size() / 2;
        int n_2 = MBmvx_2[CU32x32].size() / 2;
        MV_CT  enhanceMV_32x16_1(currCU_addx+offset32X32x[CU32x32],currCU_addy+offset32X32y[CU32x32]-8, 0, 0);
        MV_CT  enhanceMV_32x16_2(currCU_addx+offset32X32x[CU32x32],currCU_addy+offset32X32y[CU32x32]+8, 0, 0);

        nth_element(MBmvx[CU32x32].begin(), MBmvx[CU32x32].begin()+n, MBmvx[CU32x32].end()); 
        nth_element(MBmvy[CU32x32].begin(), MBmvy[CU32x32].begin()+n, MBmvy[CU32x32].end()); 
        nth_element(MBmvx_2[CU32x32].begin(), MBmvx_2[CU32x32].begin()+n_2, MBmvx_2[CU32x32].end()); 
        nth_element(MBmvy_2[CU32x32].begin(), MBmvy_2[CU32x32].begin()+n_2, MBmvy_2[CU32x32].end()); 

        if(n!=0)
        {
            enhanceMV_32x16_1.vector_X=MBmvx[CU32x32][n];
            enhanceMV_32x16_1.vector_Y=MBmvy[CU32x32][n];
            CUmv.push_back(enhanceMV_32x16_1);
        }	

        if(n_2!=0)
        {
            enhanceMV_32x16_2.vector_X=MBmvx_2[CU32x32][n_2];
            enhanceMV_32x16_2.vector_Y=MBmvy_2[CU32x32][n_2];
            CUmv.push_back(enhanceMV_32x16_2);
        }

    }//end for cu32

    //16x32
    //清空
    for(int i=0;i<4;i++)
    {
        MBmvx[i].clear();
        MBmvy[i].clear();
        MBmvx_2[i].clear();
        MBmvy_2[i].clear();
    }

    for(vector<MV_CT>::iterator it = CUmv.begin(); it != CUmv.end(); it++) 
    {
        MV_CT	getMV=*it;
        if((getMV.add_X-currCU_addx)<32&&(getMV.add_Y-currCU_addy)<32)
        { 
            if(getMV.add_X<currCU_addx+16)//比X
            {//左面
                MBmvx[0].push_back(getMV.vector_X);
                MBmvy[0].push_back(getMV.vector_Y);
            }else{//右面
                MBmvx_2[0].push_back(getMV.vector_X);
                MBmvy_2[0].push_back(getMV.vector_Y);			
            }
        }
        else if((getMV.add_X-currCU_addx)>32&&(getMV.add_Y-currCU_addy)>32)
        { 
            if(getMV.add_X<currCU_addx+48)
            {
                MBmvx[3].push_back(getMV.vector_X);
                MBmvy[3].push_back(getMV.vector_Y);
            }else{
                MBmvx_2[3].push_back(getMV.vector_X);
                MBmvy_2[3].push_back(getMV.vector_Y);			
            }
        }
        else if((getMV.add_X-currCU_addx)>32)
        { 
            if(getMV.add_X<currCU_addx+16)
            {
                MBmvx[1].push_back(getMV.vector_X);
                MBmvy[1].push_back(getMV.vector_Y);
            }else{
                MBmvx_2[1].push_back(getMV.vector_X);
                MBmvy_2[1].push_back(getMV.vector_Y);			
            }
        }
        else
        { 
            if(getMV.add_X<currCU_addx+48)
            {
                MBmvx[2].push_back(getMV.vector_X);
                MBmvy[2].push_back(getMV.vector_Y);
            }else{
                MBmvx_2[2].push_back(getMV.vector_X);
                MBmvy_2[2].push_back(getMV.vector_Y);			
            }
        }
    }


    for(int CU32x32=0;CU32x32<4;CU32x32++)
    {
        int n = MBmvx[CU32x32].size() / 2;
        int n_2 = MBmvx_2[CU32x32].size() / 2;
        MV_CT  enhanceMV_16x32_1(currCU_addx+offset32X32x[CU32x32]-8,currCU_addy+offset32X32y[CU32x32], 0, 0);
        MV_CT  enhanceMV_16x32_2(currCU_addx+offset32X32x[CU32x32]+8,currCU_addy+offset32X32y[CU32x32], 0, 0);

        nth_element(MBmvx[CU32x32].begin(), MBmvx[CU32x32].begin()+n, MBmvx[CU32x32].end()); 
        nth_element(MBmvy[CU32x32].begin(), MBmvy[CU32x32].begin()+n, MBmvy[CU32x32].end()); 
        nth_element(MBmvx_2[CU32x32].begin(), MBmvx_2[CU32x32].begin()+n_2, MBmvx_2[CU32x32].end()); 
        nth_element(MBmvy_2[CU32x32].begin(), MBmvy_2[CU32x32].begin()+n_2, MBmvy_2[CU32x32].end()); 

        if(n!=0)
        {
            enhanceMV_16x32_1.vector_X=MBmvx[CU32x32][n];
            enhanceMV_16x32_1.vector_Y=MBmvy[CU32x32][n];
            CUmv.push_back(enhanceMV_16x32_1);
        }	

        if(n_2!=0)
        {
            enhanceMV_16x32_2.vector_X=MBmvx_2[CU32x32][n_2];
            enhanceMV_16x32_2.vector_Y=MBmvy_2[CU32x32][n_2];
            CUmv.push_back(enhanceMV_16x32_2);
        }

    }//end for cu32
#endif
}


void MVmap_CT::CU32mv_mode (vector<MV_CT> &CUmv,int All_CUadd){
    int CU_offset_x = 0;
    int CU_offset_y = 0;
    int currCU_addx;
    int currCU_addy;
    int currCUadd=All_CUadd %FrameShift;
    CU_offset_x=(currCUadd)%(RowCell);
    CU_offset_y=(currCUadd-(currCUadd)%(RowCell))/RowCell;
    currCU_addx=CU_offset_x*offset64;
    currCU_addy=CU_offset_y*offset64;


    for(int CU32x32=0;CU32x32<4;CU32x32++){
        MV_CT  enhanceMV_32x32(currCU_addx+offset32X32x[CU32x32],currCU_addy+offset32X32y[CU32x32], 0, 0);
        CUmv.push_back(enhanceMV_32x32);
    }
}
