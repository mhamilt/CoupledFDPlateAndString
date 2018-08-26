%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Sanity check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% coefficients work out

Rp = k^2/(pl_rho*pl_H*pl_h^2)
Rs = k^2*st_c^2/(st_h*st_T)
Q = 2*st_T*st_h*pl_rho*pl_H*pl_h^2/(k*(st_T*st_h + pl_rho*pl_H*pl_h^2*st_c^2))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%String terms around w(1

%%% full term
tau_s*(st_lambda^2*Dn(1,:)*w1 + 2*w1(1) - 2*w2(1) + J'*(pl_mu^2*BH*u1 - 2*u1 + u2 ))

%%%current time step
tau_s*(st_lambda^2*Dn(1,:)*w1 + 2*w1(1) + J'*(pl_mu^2*BH*u1 - 2*u1)) ...
 + tau_s*(-2*w2(1) + J'*u2)

%%%previous time step
tau_s*(st_lambda^2*Dn(1,:)*w1 + 2*w1(1) - 2*w2(1) + J'*(pl_mu^2*BH*u1 - 2*u1 + u2 ))

%%String terms equivalents
tau_s*(st_lambda^2*Dn(1,:)*w1 + 2*w1(1))
tau_s*(st_lambda^2*Dn(1,:) + 2*st_I(1,:))*w1
tau_s*st_mB(1,:)*w1

%%String and plate terms equivalents current time step
tau_s*(st_lambda^2*Dn(1,:)*w1 + 2*w1(1) + J'*(pl_mu^2*BH*u1 - 2*u1))
tau_s*(st_lambda^2*Dn(1,:) + 2*st_I(1,:))*w1 + tau_s*(pl_mu^2*BH(pl_w0,:)*u1 - 2*u1(pl_w0,:))
tau_s*(st_lambda^2*Dn(1,:) + 2*st_I(1,:))*w1 + tau_s*(-pl_mB(pl_w0,:))*u1
tau_s*(st_mB(1,:)+2*st_I(1,:))*w1 + tau_s*(-pl_mB(pl_w0,:))*u1
tau_s*([st_mB(1,:)+2*st_I(1,:),-pl_mB(pl_w0,:)])*[w1;u1]

%%%NOTE
tau_s*cpl_v*[u1;w1] + tau_s*(-2*w2(1) + J'*u2)


%%%previous time step
tau_s*(-2*w2(1) + J'*u2)
tau_s*([J',-2*st_I(1,:)])*[u2;w2]
tau_s*([J',2*st_mC(1,:)])*[u2;w2]



C(ss+1,:)*[u2;w2]


%%full term
st_lambda^2*Dn(1,:)*w1 + 2*w(1) - tau_s*(st_lambda^2*Dn(1,:)*w1 + 2*w1(1) - 2*w2(1) + J'*(pl_mu^2*BH*u1 - 2*u1 + u2 ))
st_mB(1,:)*w1 + 2*w(1) - tau_s*(st_lambda^2*Dn(1,:)*w1 + 2*w1(1) - 2*w2(1) + J'*(pl_mu^2*BH*u1 - 2*u1 + u2 ))


st_lambda^2*Dn(1,:)*w1 + 2*w1(1) - tau_s*(st_lambda^2*Dn(1,:)*w1 + 2*w1(1) + J'*(pl_mu^2*BH*u1 - 2*u1))
st_mB(1,:)*w1 - tau_s*(st_lambda^2*Dn(1,:)*w1 + 2*w1(1) + (-pl_mB(pl_w0,:)*u1))
st_lambda^2*Dn(1,:)*w1 + 2*w1(1) - tau_s*cpl_v*[u1;w1]
B(ss+1,:)*[u1;w1]



st_lambda^2*Dn(1,:)*w1 + 2*w1(1) - w2(1) - tau_s*(st_lambda^2*Dn(1,:)*w1 + 2*w1(1) - 2*w2(1) + J'*(pl_mu^2*BH*u1 - 2*u1 + u2 ))
st_mB(1,:)*w1 - w2(1) - tau_s*(st_mB(1,:)*w1 - 2*w2(1) + (-pl_mB(pl_w0,:)*u1) + J'*u2 )



st_lambda^2*Dn(1,:)*w1 + 2*w1(1) - tau_s*(st_lambda^2*Dn(1,:)*w1 + 2*w1(1) + J'*(pl_mu^2*BH*u1 - 2*u1)) - w2(1) - tau_s*(-2*w2(1) + J'*u2)
B(ss+1,:)*[u1;w1] - w2(1) - tau_s*(-2*w2(1) + J'*u2)

B(ss+1,:)*[u1;w1] - C(ss+1,:)*[u2;w2]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%plate term equivalents
J'*(pl_mu^2*BH*u1 - 2*u1 + u2 )

tau_s*(st_mB(1,:)*w1-pl_mB(pl_w0,:)*u1) +  J'*u2 - 2*w2(1)


tau_s*(-pl_mB(pl_w0,:)*u1)

tau_s*([st_mB(1,:),0]*[w1;0])

tau_s*([0,-pl_mB(pl_w0,:)]*[0;u1])


tau_s*(st_mB(1,:)*w1 -pl_mB(pl_w0,:)*u1 + J'*u2)

st_lambda^2*Dn(1,:) - tau_s*(st_lambda^2*Dn(1,:) + 2*st_I(1,:) + J'*(pl_mu^2*BH*u1 - 2*pl_I(pl_w0,:)  ))
tau_s*(- 2*w2(1) + J'*u2)
% EOF
