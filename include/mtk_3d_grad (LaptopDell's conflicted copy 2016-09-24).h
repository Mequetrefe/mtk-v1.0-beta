class MTK_3DGrad {
 public:
  MTK_3DGrad();
  MTK_3DGrad(const MTK_3DGrad&);
  MTK_3DGrad(int);
  MTK_3DGrad(int, MTK_Real);
  
 private: 
    int kk_;
  MTK_Real tau_;
  MTK_DenseMatrix lap_;
  
};
