function Phi_aprimeKron=CastanedaDiazGimenezRiosRull2003_PhiaprimeMatrix(n_d,n_s,k_grid,J,zlowerbar,tauE)
% Case2_Type==2: a'(d,z,z')

N_d=prod(n_d);
N_s=prod(n_s);
k_grid=gather(k_grid);

Phi_aprimeKron=zeros([N_d,N_s,N_s]);
parfor i=1:N_d
    %Next periods assets is what you choose
    d_sub=ind2sub_homemade(n_d,i);
    d2_c=d_sub(2);
    Phi_aprimeKron_d=d2_c*ones([1,N_s,N_s]);
    %Unless you have to pay the inhertiance tax, in which case...
    if k_grid(d2_c)>=zlowerbar %if the inheritance tax is relevant
        tau_e=tauE*(k_grid(d2_c)-zlowerbar);
        aprime=k_grid(d2_c)-tau_e;
        aprime_c=dsearchn(k_grid,aprime);
        Phi_aprimeKron_d(1,J+1:2*J,1:J)=aprime_c*ones(size(Phi_aprimeKron_d(1,J+1:2*J,1:J)));
    end
    Phi_aprimeKron(i,:,:)=Phi_aprimeKron_d;
end

Phi_aprimeKron=gpuArray(Phi_aprimeKron);

% Was using the following as a check for debugging
% % Check that is seems find
% if sum(sum(sum(Phi_aprimeKron>=1)))~=N_d*N_s*N_s || sum(sum(sum(Phi_aprimeKron<=n_d(2))))~=N_d*N_s*N_s % Every choice for end-of-period assets should map to a specific tomorrows assets, for any today's shock and tomorrow's shock
%     fprintf('CastanedaDiazGimenezRiosRull2003_PhiaprimeMatrix: Seems to be a problem with Phi_aprimeKron \n')
%     sum(sum(sum(Phi_aprimeKron<1)))
%     sum(sum(sum(Phi_aprimeKron>n_d(2))))
%     min(min(min(Phi_aprimeKron)))
% end

end