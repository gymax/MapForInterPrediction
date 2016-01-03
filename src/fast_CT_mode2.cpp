#include <iostream> 
#include <fstream>
#include <vector> 
#include "fast_CT_mode2.h" 
#include <math.h>
#include <algorithm>
//block offset
/*  order       ________
    | 1  2 |
    | 3  4 |
    --------		*/
//#define  offset32X32  16
//#define  offset64X64  32
#define  offset64       64
#define  offset16       16
#define  offset64xy     32
static	int offset32x32x[] = {16,48,16,48};
static  int offset32x32y[] = {16,16,48,48};


Node::Node()
{
    is_intra=-1;
    is_skip=-1;
    partition_num=0;
    uweight=0;
    udepth=0;
    coff=-1;
    is_equalMV=-1;
}

Node::Node(int partition,int flag_intra,int flag_skip)
{
    is_intra=flag_intra;
    is_skip=flag_skip;
    partition_num=partition;
    uweight=0;
    udepth=0;
    MBadd=0;
    coff=-1;
}

//================================================================
MODEtree_CT::MODEtree_CT()
{
    dec264_tree=0;
    count_SubMB=0;
}


void MODEtree_CT::create()
{
    dec264_tree=new Node[256];//256 是基本的最大partition數
}

void MODEtree_CT::destroy()
{
    delete [] dec264_tree;
    dec264_tree=0;
}

void MODEtree_CT::show_node(int uiIdx)
{
    //cout<<<<dec264_tree[uiIdx].MBadd<<"";
    if(dec264_tree[uiIdx].MBadd!=0){
        printf("-----------MBadd= %d -----------\n",dec264_tree[uiIdx].MBadd);
        printf("partition=%d\n",dec264_tree[uiIdx].get_partnum());
        printf("weight=%d\n",dec264_tree[uiIdx].get_weight());
    }
}

void MODEtree_CT::set_node_false()
{
    Node fasle_node;
    fasle_node.set_part(-1);
    fasle_node.set_weight(-1);
    for(int i=0;i<256;i++)
        dec264_tree[i]=fasle_node;

}
//=====================================================================
MODEframe_tree_CT::MODEframe_tree_CT()
{
    dec264_CU_tree=0;

}

void MODEframe_tree_CT::create(int CU_perFrame)
{
    _CU_perFrame=CU_perFrame;
    dec264_CU_tree=new MODEtree_CT[_CU_perFrame];

    for(int i =0 ; i<_CU_perFrame; i++)
        dec264_CU_tree[i].create();


}

void MODEframe_tree_CT::destroy()
{

    for(int i =0 ; i<_CU_perFrame; i++)
        dec264_CU_tree[i].destroy();

    delete [] dec264_CU_tree;
    dec264_CU_tree=0;
}


void MODEframe_tree_CT::set_tree_dec264(std::vector<Node> dec264_node,int CUadd, int MBmap[])
{
    int uiIdx=0;
    int count_equalMV=0;
    int count_equalMV_d1[4]={0,0,0,0};
    int MB_Qtree[16]={0,1,4,5,
        2,3,6,7,
        8,9,12,13,
        10,11,14,15};
    int iMV[16][2]={0};
    int iMV_count[16]={0};
    //MB mapping for Qtree 

    for(int MBcount=0;MBcount<16;MBcount++){
        bool flag_subMB=0;
        vector <Node> ::iterator   iter; 
        iter = find_if(dec264_node.begin(),dec264_node.end(),Node_MBadd_equ(MBmap[MB_Qtree[MBcount]]));

        //
        //set depth default 16*16
        iter->set_depth(2);    
        //set weight


        if(iter->get_partnum()==1&&iter->isskip())
        {
            iter->set_depth(2); // depth 1 says skip
            iter->set_weight(0);
        }else
        {

            if(iter->get_partnum()==1)
            {
                iter->set_weight(2);
            }
            else if(iter->get_partnum()==2)
            {
                iter->set_weight(4);
            }
            else //depth ++
            {
                flag_subMB=1;
                iter->set_depth(3);
                if(!iter->isintra())
                {
                    iter->set_weight(iter->get_partnum()*64);
                }else
                {
                    iter->set_weight(-1);
                }
            }

        }
        //intra  I4很怪 不要參考
        if(iter->isintra())
        {
            iter->set_part(1);
        }

        dec264_CU_tree[CUadd].set_node(uiIdx,*iter);

        if(flag_subMB)dec264_CU_tree[CUadd].count_SubMB++;

        if(iter->isequalMV())
        {
            iMV[MBcount][0]=iter->dec_MV[0];
            iMV[MBcount][1]=iter->dec_MV[1];
        }else
        {
            iMV[MBcount][0]=65535;
            iMV[MBcount][1]=65535;	
        }


        uiIdx+=16;
    }

    //判斷是否mv相同
    //for iMV 16
    for(int i=0;i<16;i++)
    {
        if(iMV[i][0]!=65535)
        {
            for(int j=0;j<16;j++)
            {
                if(iMV[i][0]==iMV[j][0]&&iMV[i][1]==iMV[j][1])
                {
                    iMV_count[i]++;
                }

            }
        }
    }
    int MAX_iMV=0;
    for(int i=0;i<16;i++)
    {
        if(iMV_count[i]>MAX_iMV)
            MAX_iMV=iMV_count[i];
    }
    dec264_CU_tree[CUadd].equalMV=MAX_iMV;

    for(int i=0;i<16;i++)
        iMV_count[i]=0;

    //for iMV 4
    for(int k=0;k<4;k++)
    {
        for(int i=0;i<4;i++)
        {
            if(iMV[k*4+i][0]!=65535)
            {
                for(int j=0;j<4;j++)
                {
                    if(iMV[k*4+i][0]==iMV[k*4+j][0]&&iMV[k*4+i][1]==iMV[k*4+j][1])
                    {
                        iMV_count[k*4+i]++;
                    }

                }
            }
        }
        MAX_iMV=0;
        for(int i=0;i<4;i++)
        {
            if(iMV_count[k*4+i]>MAX_iMV)
                MAX_iMV=iMV_count[k*4+i];
        }
        dec264_CU_tree[CUadd].equalMV_d1[k]=MAX_iMV;
    }
    /*

       if(MBcount==0||MBcount==4||MBcount==8||MBcount==12)count_equalMV_d1=0;


       if(MBcount==0)
       {
       dec_MVx=iter->dec_MV[0];
       dec_MVy=iter->dec_MV[1];
       }

       if(iter->isequalMV())
       {
       if(dec_MVx==iter->dec_MV[0]&&dec_MVy==iter->dec_MV[1])
       {
       count_equalMV++;
       count_equalMV_d1++;
       }
       }

       if(MBcount==3||MBcount==7||MBcount==11||MBcount==15)dec264_CU_tree[CUadd].equalMV_d1[(MBcount-3)/4]=count_equalMV_d1;
       dec264_CU_tree[CUadd].equalMV=count_equalMV;
       */
}

void MODEframe_tree_CT::set_tree_dec264(int CUadd)
{
    dec264_CU_tree[CUadd].set_node_false();
}



//=====================================================================

Transcoder_CT_mode2::Transcoder_CT_mode2()
{ 
    _width=0;
    _height=0;
    _totalFrame=0; 
    total_MBA=0;
}

void Transcoder_CT_mode2::init(int width,int height,int FrameToEncode,int QP,int flag_ZB)
{ 
    _width=width;
    _height=height;
    _totalFrame=FrameToEncode;
    _flag_ZB=flag_ZB;
    _QP=QP;
    RowCell_CU=ceil(double(width/64.0));
    ColCell_CU=ceil(double(height/64.0));
    FrameShift_CU=RowCell_CU*ColCell_CU;
    _currFrame=0;
    RowCell_MB=ceil(double(width/16.0));
    ColCell_MB=ceil(double(height/16.0));
    FrameShift_MB=RowCell_MB*ColCell_MB;

    avg_MBA=0;
}

void Transcoder_CT_mode2::create()
{
    _currFrame=new MODEframe_tree_CT[_totalFrame];
    avg_MBA=new double[_totalFrame];
    for(int i =0 ; i<_totalFrame; i++)
        _currFrame[i].create(FrameShift_CU);

    //讀入黨案
    Analysis_data();

}


void Transcoder_CT_mode2::destroy()
{
    for(int i =0 ; i<_totalFrame; i++)
        _currFrame[i].destroy();

    delete [] _currFrame;
    _currFrame=0;
}



void Transcoder_CT_mode2::Analysis_data()
{

    if(_flag_ZB)
        cout<<"------------------------using ZB--------------------------\n";
    else
        cout<<"------------------------using AMBA--------------------------\n";
    ifstream fin("dec_264_mode2.txt");
    string str;
    _Count_currframe	=	0;
    int selectframe	=	0;
    int frame_end	=   0;
    double Qstep=0.625*pow(2,_QP/6.0);
    //對frame 生成_MVmap
    //讀檔

    //totalFrame-1 是因為本來就沒讀入第一張(i frame)
    for(int i = 0 ;i<_totalFrame-1 ; i++)
    {
        vector <Node> dec264_node; 
        Node curr_node;
        while(fin.good())
        {

            fin >> str;
            //判斷frame是否結束
            if((str[0]=='-')&(frame_end>0))
            {	
                avg_MBA[_Count_currframe]=total_MBA/double(FrameShift_MB);
                cout<<avg_MBA[_Count_currframe];
                total_MBA=0;
                break;  		  
            }

            //判斷是第幾個frame
            if(str[0]=='{')
            { 
                fin >> str;
                frame_end = 0;
                _Count_currframe++;
                frame_end++;
                cout << "目前是frame : " << _Count_currframe << endl;
            }
            //read  mode info


            if(str[0]=='<')
            {
                string  strMBadd;		  
                fin >> strMBadd;
                //  cout<<atoi(strMBadd.c_str())<<endl;
                curr_node.MBadd=atoi(strMBadd.c_str());

            }

            if(str[0]=='(')
            {
                string  strPartnum;
                fin >> strPartnum;
                curr_node.set_part(atoi(strPartnum.c_str()));
            }     

            if(str[0]=='!')
            {
                string  strSKIP;
                fin >> strSKIP;
                curr_node.set_skip(atoi(strSKIP.c_str()));
            }

            if(str[0]=='?')
            {
                string  strINTRA;
                fin >> strINTRA;
                curr_node.set_intra(atoi(strINTRA.c_str()));

            } 

            if(str[0]=='#')
            {
                string  str_equalMV,str_dec_MVx,str_dec_MVy;
                fin >> str_equalMV;
                curr_node.set_equalMV(atoi(str_equalMV.c_str()));
                fin >> str_dec_MVx;
                fin >> str_dec_MVy;
                curr_node.set_decMV(atoi(str_dec_MVx.c_str()),atoi(str_dec_MVy.c_str()));
            } 

            if(str[0]=='$')
            {
                string  strCOF;
                double coeff[16];
                double avg_coeff=0;
                // cout<<"------------------------using AMBA--------------------------\n";
                if(_flag_ZB)
                {
                    //  cout<<"------------------------using ZB--------------------------\n";
                    fin >> strCOF;
                    coeff[0]=atof(strCOF.c_str());
                    fin >> strCOF;
                    if(strCOF[0]=='<')
                    {
                        cout<<"讀檔錯誤"<<endl;
                        system("pause");
                    }else
                    {
                        coeff[1]=atof(strCOF.c_str());

                        for(int j=2;j<16;j++)
                        {
                            fin >> strCOF;
                            coeff[j]=atof(strCOF.c_str());
                        }

                    }


                    for(int j=0;j<16;j++)
                    {
                        if(coeff[j]<3.5*Qstep)coeff[j]=0;
                    }

                    for(int j=0;j<16;j++)
                    {
                        avg_coeff+=coeff[j];
                    }

                    curr_node.set_coff(avg_coeff);
                }
                else
                {  
                    fin >> strCOF;
                    curr_node.set_coff(atof(strCOF.c_str()));
                }

                total_MBA+=atof(strCOF.c_str()); 
                //丟進去
                dec264_node.push_back(curr_node);

                //vector dec264_node存了一個frame的mb資訊
            } 
        }//end while
        cout<<dec264_node.size()<<endl;    
        Processing_data_CU(dec264_node);



        dec264_node.clear();
    }//end for
    fin.close();

}



void Transcoder_CT_mode2::Processing_data_CU(std::vector<Node> dec264_node)
{
    int totalCU=FrameShift_CU;

    //for all CU per frame	
    for(int CUadd=0;CUadd<totalCU; CUadd++)
    {
        int MBmap[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        MBtoCU_mapping(CUadd,MBmap);

        //將mb info放入tree
        if(MBmap[0]!=-1)//-1表示邊界
        {
            _currFrame[_Count_currframe].set_tree_dec264(dec264_node,CUadd,MBmap);
        }else
        {
            _currFrame[_Count_currframe].set_tree_dec264(CUadd);
        }


    }

}

void Transcoder_CT_mode2::MBtoCU_mapping(const int CUadd, int (&MBmap)[16])
{
    //先判斷是否為邊界
    bool isBoundary =0;
    //從CUadd找到CUx,CUy
    int CU_offset_x = 0;
    int CU_offset_y = 0;
    int MB_offset_x = 0;
    int MB_offset_y = 0;
    int currCU_addx;
    int currCU_addy;
    int MBadd_LU;

    CU_offset_x=(CUadd)%(RowCell_CU);
    CU_offset_y=(CUadd-(CUadd)%(RowCell_CU))/RowCell_CU;
    currCU_addx=CU_offset_x*offset64;
    currCU_addy=CU_offset_y*offset64;

    //利用currCU_addx,currCU_addy找到MBadd_LU (Left Up)
    MB_offset_x=ceil(double(currCU_addx/16.0));
    MB_offset_y=ceil(double(currCU_addy/16.0));

    MBadd_LU=MB_offset_x+MB_offset_y*RowCell_MB;


    //+1是因為從1開始算
    //判斷右邊
    if(((CUadd+1)%RowCell_CU==0)&&_width%64!=0)
    {
        isBoundary=1;
    }

    //判斷下面
    if((CUadd+1)>(FrameShift_CU-RowCell_CU)&&_height%64!=0)
    {
        isBoundary=1;
    }

    if(!isBoundary)
    {   
        int count=0;
        //寫入MBmap
        for(int j=0;j<4;j++)
        {
            for(int i=0;i<4;i++)
            {
                MBmap[count]=MBadd_LU+i+j*RowCell_MB;
                count++;
            }
        }
    }else   //不是邊界　mapping完成
    {
        int count=0;
        //寫入MBmap
        for(int j=0;j<4;j++)
        {
            for(int i=0;i<4;i++)
            {
                MBmap[count]=-1;
                count++;
            }
        }
    //有空再做
    }


}

void Transcoder_CT_mode2::show(int frame)
{
    for(int idx=0;idx<256;idx++)
        _currFrame[frame].get_tree(12).show_node(idx);
}

MODEtree_CT Transcoder_CT_mode2::get_tree(int currframe, int currCUadd)
{
    return _currFrame[currframe].get_tree(currCUadd);
}

void Transcoder_CT_mode2::set_Tcoder_CU_AMBA(int currframe,int currCUadd)
{
    int currPIC=currframe;
    int currCU=currCUadd;
    double dec264_coff[16];
    int i=0;
    //初始化
    dec264_totalWeight=0;
    dec264_subMB=0;
    dec264_skip=0;
    dec264_AMBA_d0=0.0;

    MODEtree_CT dec264_CUtree;

    dec264_CUtree=_currFrame[currPIC].get_tree(currCU);
    for(int k=0;k<256;k+=16)//k記數
    {
        dec264_isskip[i]=dec264_CUtree.get_node(k).isskip();
        dec264_part[i]=dec264_CUtree.get_node(k).get_partnum();
        dec264_coff[i]=dec264_CUtree.get_node(k).get_coff();
        dec264_totalWeight+=dec264_CUtree.get_node(k).get_partnum();
        i++;
    }

    dec264_subMB=dec264_CUtree.count_SubMB;


    for(int k=0;k<16;k++)
    {
        dec264_AMBA_d2[k]=dec264_coff[k]/get_frameMBA(currPIC);
        if(dec264_isskip[k])
            dec264_skip++;
    }		

    for(int k=0;k<4;k++)
    {
        dec264_AMBA_d1[k]=(dec264_AMBA_d2[k*4]+dec264_AMBA_d2[k*4+1]+dec264_AMBA_d2[k*4+2]+dec264_AMBA_d2[k*4+3])/4.0;      
    }

    dec264_AMBA_d0=(dec264_AMBA_d1[0]+dec264_AMBA_d1[1]+dec264_AMBA_d1[2]+dec264_AMBA_d1[3])/4.0;

    if(dec264_part[0]==-1)//邊界
    {
        dec264_skip=-1;

    }
}


