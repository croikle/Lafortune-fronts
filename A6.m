% matrix in the 6-d wedge product space
function result = A6(eps,c,lambda,h,Z,sigma,w_star)
    A4 = A(eps,c,lambda,h,Z,sigma,w_star);
    
    function it = compute(x)
        a = A4(x);
        it = [ a(1,1) + a(2,2), a(2,3), a(2,4), -a(1,3), -a(1,4), 0 ;
               a(3,2), a(1,1) + a(3,3), a(3,4), a(1,2), 0, -a(1,4) ;
               a(4,2), a(4,3), a(1,1) + a(4,4), 0, a(1,2), a(1,3) ;
               -a(3,1), a(2,1), 0, a(2,2) + a(3,3), a(3,4), -a(2,4) ;
               -a(4,1), 0, a(2,1), a(4,3), a(2,2) + a(4,4), a(2,3) ;
               0, -a(4,1), a(3,1), -a(4,2), a(3,2), a(3,3) + a(4,4) ;
             ];
        % maybe we could generate this programmatically
    end

    result = @compute;
end
