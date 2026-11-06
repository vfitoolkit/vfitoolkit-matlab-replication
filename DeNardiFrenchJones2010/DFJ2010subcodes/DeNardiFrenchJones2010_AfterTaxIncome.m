function aftertaxincome=DeNardiFrenchJones2010_AfterTaxIncome(aprime,a,h,zeta,xi,r,earnings, taxbracket1, taxbracket2, taxbracket3, taxbracket4, taxbracket5, taxbracket6, margtaxrate0, margtaxrate1, margtaxrate2, margtaxrate3, margtaxrate4, margtaxrate5, margtaxrate6)

aftertaxincome=0;

if h==0 || h==1 % Alive (good or bad health)

    % income and taxes
    taxableincome=r*a+earnings;
    % use marginal tax rates and tax brackets
    tax=0;
    if taxableincome<taxbracket1
        tax=tax+margtaxrate0*(taxableincome-0);
    else
        tax=tax+margtaxrate0*(taxbracket1-0);
    end
    if taxableincome<taxbracket2
        tax=tax+margtaxrate1*max(taxableincome-taxbracket1,0);
    else
        tax=tax+margtaxrate1*(taxbracket2-taxbracket1);
    end
    if taxableincome<taxbracket3
        tax=tax+margtaxrate2*max(taxableincome-taxbracket2,0);
    else
        tax=tax+margtaxrate2*(taxbracket3-taxbracket2);
    end
    if taxableincome<taxbracket4
        tax=tax+margtaxrate3*max(taxableincome-taxbracket3,0);
    else
        tax=tax+margtaxrate3*(taxbracket4-taxbracket3);
    end
    if taxableincome<taxbracket5
        tax=tax+margtaxrate4*max(taxableincome-taxbracket4,0);
    else
        tax=tax+margtaxrate4*(taxbracket5-taxbracket4);
    end
    if taxableincome<taxbracket6
        tax=tax+margtaxrate5*max(taxableincome-taxbracket5,0);
    else
        tax=tax+margtaxrate5*(taxbracket6-taxbracket5);
    end
    if taxableincome>=taxbracket6 % top bracket
        tax=tax+margtaxrate6*max(taxableincome-taxbracket6,0);
    end
    % after-tax income is then just
    aftertaxincome=taxableincome-tax;

end


end