

class MTK_3DDiv {
 public:
  MTK_3DDiv();
  MTK_3DDiv(const MTK_3DDiv&);
  MTK_3DDiv(int);
  MTK_3DDiv(int, MTK_Real);
  MTK_3DDiv(int, MTK_Real,int num_xx,int num_yy);
  
 private: 
    int kk_;
  MTK_Real tau_;
  MTK_DenseMatrix lap_;
  
};