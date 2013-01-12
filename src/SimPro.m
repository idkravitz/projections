

function [z, reps, iter, lmb, kvec, R1, info] = SimPro(X, maxit, calceps, verbose, kvec0, R)

if verbose(1) >= 0
        printf("\n%s\n\n", "$Id: SimPro.w,v 2.20 2011/11/17 12:40:14 nurmi Exp nurmi $");
        printf(" dim = %4d nvec = %5d\n", size(X));
endif


OPTIMAL                 =   0;
OVERSIZED_BASIS         = -19;
LINDEP_BASIS            = -20;
INIT_ERROR              = -21;
NO_WAY_TOI              = -22;
CHOLIN_CRASH            = -23;
BAD_OLD_L               = -24;
MULTI_IN                = -25;
DBLE_VTX                = -26;
EXCESS_ITER             = -27;
RUNNING                 = -28;
NEGGRAMM                = -29;
info = RUNNING;



macheps = 1.0e-16; 

reps = 0;
iter = [ 0 0 ];
[ dim nvec ] = size(X);
if size(kvec0) == 0
        
        for j = 1:nvec
                gg(j) = dot(X(:,j), X(:,j));
        endfor
        [ vx ix ] = min(gg);
        kvec = [ix]; lmb = lmb_old = [1]; baslen = 1; z = X(:,ix);
        R1 = R = sqrt(vx)*ones(1,1);
        if vx < calceps*calceps
                reps = vx;
                info = OPTIMAL;
                return
        endif
        
else
        
        R1 = R;
        kvec = kvec0;
        lmb = lmb_old = bari(X, kvec, R);
        baslen = columns(kvec);

        bad_old = lmb_old < 0;
        if any(bad_old)
                idx_bad_old = [1:rows(lmb_old)](bad_old)
                printf(" XXX WOW! Bad old lambda !\n");
                for i = idx_bad_old
                        printf(" i %4d lmb_old %12.6e\n", i, lmb_old(i));
                endfor
                info = BAD_OLD_L;
                return;
        endif
        
        if !all(lmb_old + calceps/baslen > 0)
                printf(" Initialization error. Negative baricentrics.\n");
                info = INIT_ERROR;
                return
        endif
        [ lmb_old R ins ] = cleanbas(lmb_old, R, 1.e-08);
        lmb_old /= sum(lmb_old);

        bad_old = lmb_old < 0;
        if any(bad_old)
                idx_bad_old = [1:rows(lmb_old)](bad_old)
                printf(" XXX WOW! Bad old lambda !\n");
                for i = idx_bad_old
                        printf(" i %4d lmb_old %12.6e\n", i, lmb_old(i));
                endfor
                info = BAD_OLD_L;
                return;
        endif
        
        kvec = kvec(ins);
        baslen = columns(kvec);
        z = X(:, kvec)*lmb_old;
        R1 = R;
        
endif
info = RUNNING;


while ( iter(1) < maxit )
        


        zz = z'*z;
        dkvec = setdiff([1:nvec], abs(kvec));
        vcos = z'*X(:, dkvec) - zz;
        if all(vcos + calceps > 0)
                retcode = OPTIMAL;
                reps = min(vcos);
                kvec = abs(kvec);
                info = OPTIMAL;
                break;
        endif
        

        [reps, iv] = min(vcos);
        if columns(iv) > 1
                printf(" Multiple min-vec. I am not ready to handle this case.\n");
                info = MULTI_IN;
                return;
        endif
        ix = dkvec(iv);
        kvec = [ kvec ix ];
        baslen++;
        

        xi = X(:, ix);
        [ R1 cinfo ] = cholinsert(R, baslen, X(:, abs(kvec))'*xi);
        switch (cinfo)
        case 0
                
                lmb_new = R1 \ ( R1' \ ones(baslen, 1));
                lmb_new /= sum(lmb_new);
                R = R1;
                
        case 1
                kvec(baslen) *= -1;
                lmb_new = bari(X, kvec, R);
                
                R1 = R;
        case 2
                printf(" Smth is TERRIBLY wrong.\n");
                info = CHOLIN_CRASH;
                return;
        endswitch
        

        verbose(1) > 0 && mod(iter(1), verbose(1)) == 0 &&                      \
                printf(" +++ iter %4d in %4d basis %4d zx %12.4e zz %16.8e\n",  \
                   iter(1), ix, baslen, reps, zz);
        
        
        first = 1;
        while 1
                

                if all(lmb_new + calceps/rows(lmb_new) > 0)
                        z = X(:, abs(kvec))*lmb_new;
                        lmb_old = lmb_new;
                        
                        bad_old = lmb_old < 0;
                        if any(bad_old)
                                idx_bad_old = [1:rows(lmb_old)](bad_old)
                                printf(" XXX WOW! Bad old lambda !\n");
                                for i = idx_bad_old
                                        printf(" i %4d lmb_old %12.6e\n", i, lmb_old(i));
                                endfor
                                info = BAD_OLD_L;
                                return;
                        endif
                        
                        break
                endif
                

                if rows(lmb_old) < rows(lmb_new)
                        lmb_old = [ lmb_old; 0 ];
                        if rows(lmb_old) < rows(lmb_new)
                                printf("Somehow we jump 2 dimensions up !\n");
                                info = DBLE_VTX;
                                return
                        endif
                endif

                zl = lmb_old - lmb_new;
                imin = (zl > 0) & (lmb_old > 0) & (lmb_new < 0);
                if imin == 0
                        printf("No way to decrease z. Pass away ...\n");
                        printf("POST MORTUM:\n");
                        zl_pos = [1:rows(lmb_old)](zl > 0)
                        lmb_old_neg = [1:rows(lmb_old)](lmb_old < 0);
                        disp([ lmb_old_neg lmb_old(lmb_old_neg) ]);
                        
                        bad_old = lmb_old < 0;
                        if any(bad_old)
                                idx_bad_old = [1:rows(lmb_old)](bad_old)
                                printf(" XXX WOW! Bad old lambda !\n");
                                for i = idx_bad_old
                                        printf(" i %4d lmb_old %12.6e\n", i, lmb_old(i));
                                endfor
                                info = BAD_OLD_L;
                                return;
                        endif
                        
                        lmb_new_neg = [1:rows(lmb_old)](lmb_new < -calceps);
                        badidx = (lmb_old < calceps) & (lmb_new > -calceps);
                        iox = [1: rows(lmb_old)](badidx)
                        lmb = lmb_old;
                        info = NO_WAY_TOI;
                        return;
                endif
                [t, inx] = min( lmb_old(imin) ./ zl(imin) );
                
                 
                lmb_t = t*lmb_new + (1-t)*lmb_old;
                iox = [1:baslen](imin)(inx);
                kout = kvec(iox);
                kvec_new = [ kvec(1 : (iox-1)) kvec((iox+1) : columns(kvec) ) ];
                lmb_t = [ lmb_t'(1 : (iox-1)) lmb_t'((iox+1) : rows(lmb_t) ) ]';
                kvec = kvec_new; baslen--;
                zt = X(:, abs(kvec))*lmb_t;
                

                R1 = choldelete (R, iox);
                R = R1;
                if kvec(baslen) < 0
                        kvec(baslen) *= -1;
                        [ R1 infor ] = cholinsert(R, baslen, X(:, abs(kvec))'*X(:, kvec(baslen)) );
                        if infor == 0
                                R = R1;
                        else
                                
                                printf(" XXX Gramm matrix is not positive definite.");
                                [ vmin mx ] = min(diag(R1));
                                printf(" min %12.4e mx %5d. Refactor.\n", vmin, mx);
                                t0 = time();
                                [ R pr ] = chol(X(:, kvec)'*X(:, kvec));
                                printf(" Time to refactor %8.2f sec.\n", time()-t0);
                                if pr > 0
                                        printf(" XXX Refactoring failed. Continue at your own risk.\n");
                                        info = NEGGRAMM;
                                endif
                                
                        endif
                endif
               % [ lmb_t R ins ] = cleanbas(lmb_t, R, 1.e-08);
%ins
%size(R)
%lmb_t
               % kvec = kvec(ins);
                R1 = R;
                baslen = columns(kvec);
                
                

                lmb_new = bari(X, kvec, R);
                lmb_old = lmb_t;

                bad_old = lmb_old < 0;
                if any(bad_old)
                        idx_bad_old = [1:rows(lmb_old)](bad_old)
                        printf(" XXX WOW! Bad old lambda !\n");
                        for i = idx_bad_old
                                printf(" i %4d lmb_old %12.6e\n", i, lmb_old(i));
                        endfor
                        info = BAD_OLD_L;
                        return;
                endif
                
                
                iter(2)++;

                verbose(2) > 0 && mod(iter(2), verbose(2)) == 0 &&                      \
                        printf(" --- intc %4d out  %4d (%4d) basis %4d step %12.4e zz %14.6e\n", \
                           iter(2), kout, iox, baslen, t, zz);
                
                t <= 0 && verbose(1) >= 0 && printf("XXX step t negative %e\n", t);
                
                first = 0;
        endwhile
        
        iter(1)++;
endwhile
(iter(1) >= maxit) && (info = EXCESS_ITER);


if(verbose(1) >= 0)
        reps = min(z'*X - z'*z);
        printf("\n-> Status               %6d\n", info);
        printf(  "-> External iterations  %6d\n", iter(1));
        printf(  "-> Internal iterations  %6d\n", iter(2));
        printf(  "-> Basis size %6d\n", columns(kvec));
        printf(  "-> Optimality %15.4e\n", abs(reps));
        printf(  "-> Distance   %20.8e\n\n", sqrt(zz));
endif


bad_old = lmb_old < 0;
if any(bad_old)
        idx_bad_old = [1:rows(lmb_old)](bad_old)
        printf(" XXX WOW! Bad old lambda !\n");
        for i = idx_bad_old
                printf(" i %4d lmb_old %12.6e\n", i, lmb_old(i));
        endfor
        info = BAD_OLD_L;
        return;
endif

lmb = lmb_old;

endfunction


function [ lmb_clean RC isn ] = cleanbas(lmb, R, calceps)
RC = R; lmb_clean = lmb; isn = [1 : rows(lmb)];
is = [1 : rows(lmb)](abs(lmb) < calceps/rows(lmb));
if columns(is) == 0
        return
is
isn = setdiff(1:rows(lmb), is);
lmb_clean = lmb(isn);
endif
for i = columns(is) - [ 1 : columns(is)] + 1
        RC = choldelete(RC, is(i));
endfor
endfunction


function [ lmb ] = bari(X, kvec, R)
baslen = columns(kvec);
if kvec(baslen) > 0

        lmb = R \ ( R' \ ones(baslen, 1));
        lmb /= sum(lmb);
        
        return;
endif
 
kin = abs(kvec(baslen)); kv = kvec(1:baslen-1);
rhs = (X(:,kv))'*X(:, kin);
z = R \ ( R' \ rhs );
q = sum(z);
if abs(1-q) < 1.e-10 
   printf(" Problem in full basis (newbari).");
   printf(" Continue at your own risk.\n");
endif
q = 1/(1 - q);
lmb = q*[ -z; 1 ];
lmb /= sum(lmb);
 
endfunction
