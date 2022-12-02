# Global to Local: A Hierarchical Detection Algorithm for Hyperspectral Image Target Detection
----------
## Zhonghao Chen, Zhengtao Lu, Hongmin Gao, Yiyan Zhang, Danfeng Hong, Bing Zhang
----------
The code in this toolbox implements the "Global to Local: A Hierarchical Detection Algorithm for Hyperspectral Image Target Detection". More specifically, it is detailed as follow.
## System-specific notes
----------
The data were generated by Matlab R2016b or higher versions. 

## How to use it?
----------

This toolbox consists of two modules, i.e, extended morphological attribute profile (EMAP) and diverse-direction constrained energy minimization (D2CEM) detector. EMAP is run by **EMAP.m** in the ‘./EMAP’ folder, and its support functions are in the ‘./EMAP/supportFunc’ folder. D2CEM is run by **D2CEM.m** in the ‘./G2LHTD’ folder, and the **func_pear.m**, **func_calCorrMat.m** and **hyperSam.m** are used to support the calculation of D2CEM detector.
Here an example experiment is given by using **San Diego hyperspectral data**. Directly run **detect.m** in the ‘./G2LHTD’ folder with different parameter settings to produce the results. For more details, please  do not hesitate to contact us.


## Citation
----------
Please kindly cite the papers if this code is useful and helpful for your research.\
@ARTICLE{9968036,\
  author={Chen, Zhonghao and Lu, Zhengtao and Gao, Hongmin and Zhang, Yiyan and Zhao, Jia and Hong, Danfeng and Zhang, Bing},\
  journal={IEEE Transactions on Geoscience and Remote Sensing}, \
  title={Global to Local: A Hierarchical Detection Algorithm for Hyperspectral Image Target Detection}, \
  year={2022}, \
  volume={}, \
  number={}, \ 
  pages={1-1},\ 
  doi={10.1109/TGRS.2022.3225902}}
  
  ## Contact Information
  ----------
  Zhonghao Chen: chenzhonghao@hhu.edu.cn\
  Zhengtao Lu: luzhengtao@hhu.edu.cn
  
