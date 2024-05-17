function [Bus_Line] = z_matrix_Ieee13(n_branch)


A12 = [   0.0076 + 0.0223i   0.0034 + 0.0110i   0.0035 + 0.0093i
   0.0034 + 0.0110i   0.0074 + 0.0229i   0.0034 + 0.0084i
   0.0035 + 0.0093i   0.0034 + 0.0084i   0.0075 + 0.0227i];
%%
A23 = [ 0  0                  0
        0  0.0724 + 0.0743i   0.0113 + 0.0251i
        0  0.0113 + 0.0251i   0.0727 + 0.0737i];
%%
A24 = [   0.0253 + 0.0743i   0.0114 + 0.0366i   0.0115 + 0.0309i
   0.0114 + 0.0366i   0.0246 + 0.0765i   0.0112 + 0.0281i
   0.0115 + 0.0309i   0.0112 + 0.0281i   0.0249 + 0.0755i];
%%
A25  = [   0.0412 + 0.0646i   0.0086 + 0.0232i   0.0085 + 0.0275i
   0.0086 + 0.0232i   0.0409 + 0.0656i   0.0084 + 0.0211i
   0.0085 + 0.0275i   0.0084 + 0.0211i   0.0407 + 0.0663i];
%%
A46  = [   0.0506 + 0.1485i   0.0228 + 0.0732i   0.0231 + 0.0618i
   0.0228 + 0.0732i   0.0492 + 0.1529i   0.0224 + 0.0562i
   0.0231 + 0.0618i   0.0224 + 0.0562i   0.0498 + 0.1510i];
%%
A57 = [   0.0303 + 0.0891i   0.0137 + 0.0439i   0.0138 + 0.0371i
   0.0137 + 0.0439i   0.0295 + 0.0917i   0.0134 + 0.0337i
   0.0138 + 0.0371i   0.0134 + 0.0337i   0.0299 + 0.0906i];
%%
A38 = [0   0                  0
       0   0.0435 + 0.0446i   0.0068 + 0.0151i
       0   0.0068 + 0.0151i   0.0436 + 0.0442i];
%%
A69 = 1.0e-03 *[   0.0758 + 0.2228i   0.0341 + 0.1098i   0.0346 + 0.0927i
   0.0341 + 0.1098i   0.0739 + 0.2294i   0.0336 + 0.0842i
   0.0346 + 0.0927i   0.0336 + 0.0842i   0.0747 + 0.2265i];
%%
A610 = [   0.0379 + 0.1114i   0.0171 + 0.0549i   0.0173 + 0.0464i
   0.0171 + 0.0549i   0.0369 + 0.1147i   0.0168 + 0.0421i
   0.0173 + 0.0464i   0.0168 + 0.0421i   0.0374 + 0.1133i];
%%
A611 = [   0.0435 + 0.0446i  0   0.0068 + 0.0151i
                          0  0   0
           0.0068 + 0.0151i  0    0.0436 + 0.0442i];
%%
A912 = [   0.0433 + 0.0240i   0.0174 + 0.0015i   0.0155 - 0.0010i
   0.0174 + 0.0015i   0.0428 + 0.0217i   0.0174 + 0.0015i
   0.0155 - 0.0010i   0.0174 + 0.0015i   0.0433 + 0.0240i];
%%
A1113 = [0  0  0
         0  0  0
         0  0  0.0436 + 0.0442i];
%%
A1114 = [ 0.1175 + 0.0449i     0     0
                         0     0     0
                         0     0     0];

%%

Bus_Line = cell(1,n_branch);

Bus_Line(1,1) = {A12};    Bus_Line(1,2) = {A23};    Bus_Line(1,3) = {A24};
Bus_Line(1,4) = {A25};    Bus_Line(1,5) = {A46};    Bus_Line(1,6) = {A57};
Bus_Line(1,7) = {A38};    Bus_Line(1,8) = {A69};    Bus_Line(1,9) = {A610};
Bus_Line(1,10) = {A611};  Bus_Line(1,11) = {A912};  Bus_Line(1,12) = {A1113};
Bus_Line(1,13) = {A1114};

end