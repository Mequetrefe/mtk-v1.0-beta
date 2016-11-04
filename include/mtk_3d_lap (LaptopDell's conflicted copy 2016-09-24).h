class MTK_3DLap {
 public:
  MTK_3DLap();
  MTK_3DLap(const MTK_3DLap&);
  MTK_3DLap(int);
  MTK_3DLap(int, MTK_Real);
  
 private: 
  int kk_;
  MTK_Real tau_;
  MTK_DenseMatrix lap_;
};
