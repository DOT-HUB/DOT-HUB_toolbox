function [refpts_10_5, refpts_10_5_txt] = DOTHUB_getTen5points(surf,refpts)
%
% USAGE:
%
%    [refpts_10_5 refpts_10_5_txt] = get_10_5_pts(surf,refpts)
%
% INPUTS:
%    
%    surf    - Array of vertices of the head surface. 
%
%    refpts  - A set of 5 landmark points in the order: Nz, Iz, Ar, Al, Cz
%              where Cz is an initial guess (here referred to as Czi) 
%              marking some position close to the top of the head surface
%              passed as the first argument.
%
% OUTPUTS:
%     
%    refpts_10_5     - ...x3 matrix of coordinates of the 10-5 points.
%                       (Doesn't include initial reference points)
%
%    refpts_10_5_txt - Names each anatomical point in refpts_10_5 
%                       according to the 10-5 convention. 
%
% DESCRIPTION:
%
%    User picks the reference points Nz, Iz, Ar, Al and Cz (Cz can be less exact 
%    than the others, basically anywhere near the top/center of the head),
%    on a scan of the head.
%
%    *** NOTE *** this algorithm needs two conditions to work properly: 
%    a) the mesh density has to be reasonably high (a bit of trial and error),  
%    and b) the reference points should be on actual head voxels rather than 
%    air (or some other medium). 
%
% 
% EXAMPLE:
%
%    hseg = uint8(load_image('mri/hseg.bin'));
%    [f v] = isosurface(hseg,0.9);
%    v = [v(:,2) v(:,1) v(:,3)];
%    load refpts.txt -ascii
%    [refpts_10_5 refpts_10_5_txt] = get_10_5_pts(v, refpts);
%
%
% AUTHOR: Jay Dubb (jdubb@nmr.mgh.harvard.edu)
% DATE:   05/07/2009
%
% Modified by: S. Brigadoi 24/04/2013: compute the 10-5 positions instead
%                                      of the 10-20
%                                      and the procedure to choose the best
%                                      curve where to set the 10-5 points
%                                      has been changed. 
%
 
    Nz = refpts(1,:);
    Iz = refpts(2,:);
    Ar = refpts(3,:);
    Al = refpts(4,:);
    Czi = refpts(5,:);

    % Distance threshold for the plane_equation to determine 
    % the surface points on the head which intersect with a plane
    dt=.6;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find intersection of Al-Czi-Ar plane and mesh surface
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find the midpoint Czi of curve from Al to Ar
    plane = plane_equation(Al, Ar, Czi);
    [curve_pts_AlCziAr,len_AlCziAr] = curveGen(Al, Ar, Czi, plane, surf, dt);
    %%display(['Initial Al-Ar curve length: ' num2str(len_AlCziAr)]);
    
    %%display('Making initial guess at Cz point.');
    % Recalculate Czi - we narrow down the location of Cz
    Czi = curveWalk(curve_pts_AlCziAr, len_AlCziAr/2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find intersection of Nz-Czi-Iz plane and mesh surface
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plane = plane_equation(Nz, Iz, Czi);
    [curve_pts_NzCziIz,len_NzCziIz] = curveGen(Nz, Iz, Czi, plane, surf, dt);
    %%display(['Nz-Czi-Iz curve length: ' num2str(len_NzCziIz)]);
    
    [Cz,~] = curveWalk(curve_pts_NzCziIz, len_NzCziIz/2);
    % Every point is then brought back to the surface because it is
    % computed on the interpolated curve that has points that do not
    % belong to the surface. Thus in order to have the 10-5 position on the real
    % mesh (its node), every point computed is forced to the surface mesh
    Cz = bringPtsToSurf(surf,Cz); 
    %%display(['Cz = [' num2str(Cz(1,1)) ' ' num2str(Cz(1,2)) ' ' num2str(Cz(1,3)) ']']);
    
    [NFpz,~] = curveWalk(curve_pts_NzCziIz, 0.5*len_NzCziIz/10);
    NFpz = bringPtsToSurf(surf,NFpz);
    %%display(['NFpz = [' num2str(NFpz(1,1)) ' ' num2str(NFpz(1,2)) ' ' num2str(NFpz(1,3)) '], ' num2str(len) ' away from Nz']);

    [Fpz,~] = curveWalk(curve_pts_NzCziIz, 1*len_NzCziIz/10);
    Fpz = bringPtsToSurf(surf,Fpz);
    %%display(['Fpz = [' num2str(Fpz(1,1)) ' ' num2str(Fpz(1,2)) ' ' num2str(Fpz(1,3)) '], ' num2str(len) ' away from Nz']);

    [AFpz,~] = curveWalk(curve_pts_NzCziIz, 1.5*len_NzCziIz/10);
    AFpz = bringPtsToSurf(surf,AFpz);
    %%display(['AFpz = [' num2str(AFpz(1,1)) ' ' num2str(AFpz(1,2)) ' ' num2str(AFpz(1,3)) '], ' num2str(len) ' away from Nz']);
    
    [AFz,~] = curveWalk(curve_pts_NzCziIz, 2*len_NzCziIz/10);
    AFz = bringPtsToSurf(surf,AFz);
    %%display(['AFz = [' num2str(AFz(1,1)) ' ' num2str(AFz(1,2)) ' ' num2str(AFz(1,3)) '], ' num2str(len) ' away from Nz']);

    [AFFz,~] = curveWalk(curve_pts_NzCziIz, 2.5*len_NzCziIz/10);
    AFFz = bringPtsToSurf(surf,AFFz);
    %display(['AFFz = [' num2str(AFFz(1,1)) ' ' num2str(AFFz(1,2)) ' ' num2str(AFFz(1,3)) '], ' num2str(len) ' away from Nz']);
    
    [Fz,~]  = curveWalk(curve_pts_NzCziIz, 3*len_NzCziIz/10);
    Fz = bringPtsToSurf(surf,Fz);
    %display(['Fz = [' num2str(Fz(1,1)) ' ' num2str(Fz(1,2)) ' ' num2str(Fz(1,3)) '], ' num2str(len) ' away from Nz']);

    [FFCz,~] = curveWalk(curve_pts_NzCziIz, 3.5*len_NzCziIz/10);
    FFCz = bringPtsToSurf(surf,FFCz);
    %display(['FFCz = [' num2str(FFCz(1,1)) ' ' num2str(FFCz(1,2)) ' ' num2str(FFCz(1,3)) '], ' num2str(len) ' away from Nz']);
    
    [FCz,~] = curveWalk(curve_pts_NzCziIz, 4*len_NzCziIz/10);
    FCz = bringPtsToSurf(surf,FCz);
    %display(['FCz = [' num2str(FCz(1,1)) ' ' num2str(FCz(1,2)) ' ' num2str(FCz(1,3)) '], ' num2str(len) ' away from Nz']);
    
    [FCCz,~] = curveWalk(curve_pts_NzCziIz, 4.5*len_NzCziIz/10);
    FCCz = bringPtsToSurf(surf,FCCz);
    %display(['FCCz = [' num2str(FCCz(1,1)) ' ' num2str(FCCz(1,2)) ' ' num2str(FCCz(1,3)) '], ' num2str(len) ' away from Nz']);

    [CCPz,~] = curveWalk(curve_pts_NzCziIz, 5.5*len_NzCziIz/10);
    CCPz = bringPtsToSurf(surf,CCPz);
    %display(['CCPz = [' num2str(CCPz(1,1)) ' ' num2str(CCPz(1,2)) ' ' num2str(CCPz(1,3)) '], ' num2str(len) ' away from Nz']);
    
    [CPz,~] = curveWalk(curve_pts_NzCziIz, 6*len_NzCziIz/10);
    CPz = bringPtsToSurf(surf,CPz);
    %display(['CPz = [' num2str(CPz(1,1)) ' ' num2str(CPz(1,2)) ' ' num2str(CPz(1,3)) '], ' num2str(len) ' away from Nz']);
    
    [CPPz,~] = curveWalk(curve_pts_NzCziIz, 6.5*len_NzCziIz/10);
    CPPz = bringPtsToSurf(surf,CPPz);
    %display(['CPPz = [' num2str(CPPz(1,1)) ' ' num2str(CPPz(1,2)) ' ' num2str(CPPz(1,3)) '], ' num2str(len) ' away from Nz']);

    [Pz,~]  = curveWalk(curve_pts_NzCziIz, 7*len_NzCziIz/10);
    Pz = bringPtsToSurf(surf,Pz);
    %display(['Pz = [' num2str(Pz(1,1)) ' ' num2str(Pz(1,2)) ' ' num2str(Pz(1,3)) '], ' num2str(len) ' away from Nz'])

    [PPOz,~] = curveWalk(curve_pts_NzCziIz, 7.5*len_NzCziIz/10);
    PPOz = bringPtsToSurf(surf,PPOz);
    %display(['PPOz = [' num2str(PPOz(1,1)) ' ' num2str(PPOz(1,2)) ' ' num2str(PPOz(1,3)) '], ' num2str(len) ' away from Nz']);
    
    [POz,~] = curveWalk(curve_pts_NzCziIz, 8*len_NzCziIz/10);
    POz = bringPtsToSurf(surf,POz);
    %display(['POz = [' num2str(POz(1,1)) ' ' num2str(POz(1,2)) ' ' num2str(POz(1,3)) '], ' num2str(len) ' away from Nz']);
    
    [POOz,~] = curveWalk(curve_pts_NzCziIz, 8.5*len_NzCziIz/10);
    POOz = bringPtsToSurf(surf,POOz);
    %display(['POOz = [' num2str(POOz(1,1)) ' ' num2str(POOz(1,2)) ' ' num2str(POOz(1,3)) '], ' num2str(len) ' away from Nz']);

    [Oz,~]  = curveWalk(curve_pts_NzCziIz, 9*len_NzCziIz/10);
    Oz = bringPtsToSurf(surf,Oz);
    %display(['Oz = [' num2str(Oz(1,1)) ' ' num2str(Oz(1,2)) ' ' num2str(Oz(1,3)) '], ' num2str(len) ' away from Nz']);
    
    [OIz,~]  = curveWalk(curve_pts_NzCziIz, 9.5*len_NzCziIz/10);
    OIz = bringPtsToSurf(surf,OIz);
    %display(['OIz = [' num2str(OIz(1,1)) ' ' num2str(OIz(1,2)) ' ' num2str(OIz(1,3)) '], ' num2str(len) ' away from Nz']);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find point along Al-Cz-Ar curve 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Recompute curve Al to Ar and its length
    plane = plane_equation(Al, Ar, Cz);
    [curve_pts_AlCzAr,len_AlCzAr] = curveGen(Al, Ar, Cz, plane, surf, dt);
    %display(['Al-Cz-Ar curve length: ' num2str(len_AlCzAr)]);

    [T9h,~] = curveWalk(curve_pts_AlCzAr, 0.5*len_AlCzAr/10);
    T9h = bringPtsToSurf(surf,T9h);
    %display(['T9h = [' num2str(T9h(1,1)) ' ' num2str(T9h(1,2)) ' ' num2str(T9h(1,3)) '], ' num2str(len) ' away from Al']);
    
    [T7,~] = curveWalk(curve_pts_AlCzAr, 1*len_AlCzAr/10);
    T7 = bringPtsToSurf(surf,T7);
    %display(['T7 = [' num2str(T7(1,1)) ' ' num2str(T7(1,2)) ' ' num2str(T7(1,3)) '], ' num2str(len) ' away from Al']);
    
    [T7h,~] = curveWalk(curve_pts_AlCzAr, 1.5*len_AlCzAr/10);
    T7h = bringPtsToSurf(surf,T7h);
    %display(['T7h = [' num2str(T7h(1,1)) ' ' num2str(T7h(1,2)) ' ' num2str(T7h(1,3)) '], ' num2str(len) ' away from Al']);

    [C5,~] = curveWalk(curve_pts_AlCzAr, 2*len_AlCzAr/10);
    C5 = bringPtsToSurf(surf,C5);
    %display(['C5 = [' num2str(C5(1,1)) ' ' num2str(C5(1,2)) ' ' num2str(C5(1,3)) '], ' num2str(len) ' away from Al']);
    
    [C5h,~] = curveWalk(curve_pts_AlCzAr, 2.5*len_AlCzAr/10);
    C5h = bringPtsToSurf(surf,C5h);
    %display(['C5h = [' num2str(C5h(1,1)) ' ' num2str(C5h(1,2)) ' ' num2str(C5h(1,3)) '], ' num2str(len) ' away from Al']);

    [C3,~] = curveWalk(curve_pts_AlCzAr, 3*len_AlCzAr/10);
    C3 = bringPtsToSurf(surf,C3);
    %display(['C3 = [' num2str(C3(1,1)) ' ' num2str(C3(1,2)) ' ' num2str(C3(1,3)) '], ' num2str(len) ' away from Al']);
    
    [C3h,~] = curveWalk(curve_pts_AlCzAr, 3.5*len_AlCzAr/10);
    C3h = bringPtsToSurf(surf,C3h);
    %display(['C3h = [' num2str(C3h(1,1)) ' ' num2str(C3h(1,2)) ' ' num2str(C3h(1,3)) '], ' num2str(len) ' away from Al']);

    [C1,~] = curveWalk(curve_pts_AlCzAr, 4*len_AlCzAr/10);
    C1 = bringPtsToSurf(surf,C1);
    %display(['C1 = [' num2str(C1(1,1)) ' ' num2str(C1(1,2)) ' ' num2str(C1(1,3)) '], ' num2str(len) ' away from Al']);
    
    [C1h,~] = curveWalk(curve_pts_AlCzAr, 4.5*len_AlCzAr/10);
    C1h = bringPtsToSurf(surf,C1h);
    %display(['C1h = [' num2str(C1h(1,1)) ' ' num2str(C1h(1,2)) ' ' num2str(C1h(1,3)) '], ' num2str(len) ' away from Al']);

    [C2h,~] = curveWalk(curve_pts_AlCzAr, 5.5*len_AlCzAr/10);
    C2h = bringPtsToSurf(surf,C2h);
    %display(['C2h = [' num2str(C2h(1,1)) ' ' num2str(C2h(1,2)) ' ' num2str(C2h(1,3)) '], ' num2str(len) ' away from Al']);
    
    [C2,~] = curveWalk(curve_pts_AlCzAr, 6*len_AlCzAr/10);
    C2 = bringPtsToSurf(surf,C2);
    %display(['C2 = [' num2str(C2(1,1)) ' ' num2str(C2(1,2)) ' ' num2str(C2(1,3)) '], ' num2str(len) ' away from Al']);

    [C4h,~] = curveWalk(curve_pts_AlCzAr, 6.5*len_AlCzAr/10);
    C4h = bringPtsToSurf(surf,C4h);
    %display(['C4h = [' num2str(C4h(1,1)) ' ' num2str(C4h(1,2)) ' ' num2str(C4h(1,3)) '], ' num2str(len) ' away from Al']);
    
    [C4,~] = curveWalk(curve_pts_AlCzAr, 7*len_AlCzAr/10);
    C4 = bringPtsToSurf(surf,C4);
    %display(['C4 = [' num2str(C4(1,1)) ' ' num2str(C4(1,2)) ' ' num2str(C4(1,3)) '], ' num2str(len) ' away from Al']);

    [C6h,~] = curveWalk(curve_pts_AlCzAr, 7.5*len_AlCzAr/10);
    C6h = bringPtsToSurf(surf,C6h);
    %display(['C6h = [' num2str(C6h(1,1)) ' ' num2str(C6h(1,2)) ' ' num2str(C6h(1,3)) '], ' num2str(len) ' away from Al']);
    
    [C6,~] = curveWalk(curve_pts_AlCzAr, 8*len_AlCzAr/10);
    C6 = bringPtsToSurf(surf,C6);
    %display(['C6 = [' num2str(C6(1,1)) ' ' num2str(C6(1,2)) ' ' num2str(C6(1,3)) '], ' num2str(len) ' away from Al']);

    [T8h,~] = curveWalk(curve_pts_AlCzAr, 8.5*len_AlCzAr/10);
    T8h = bringPtsToSurf(surf,T8h);
    %display(['T8h = [' num2str(T8h(1,1)) ' ' num2str(T8h(1,2)) ' ' num2str(T8h(1,3)) '], ' num2str(len) ' away from Al']);
    
    [T8,~] = curveWalk(curve_pts_AlCzAr, 9*len_AlCzAr/10);
    T8 = bringPtsToSurf(surf,T8);
    %display(['T8 = [' num2str(T8(1,1)) ' ' num2str(T8(1,2)) ' ' num2str(T8(1,3)) '], ' num2str(len) ' away from Al']);
    
    [T10h,~] = curveWalk(curve_pts_AlCzAr, 9.5*len_AlCzAr/10);
    T10h = bringPtsToSurf(surf,T10h);
    %display(['T10h = [' num2str(T10h(1,1)) ' ' num2str(T10h(1,2)) ' ' num2str(T10h(1,3)) '], ' num2str(len) ' away from Al']);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find point along Fpz-T7-Oz curve 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    plane = plane_equation(Fpz, T7, Oz);
    [curve_pts_FpzT7Oz,len_FpzT7Oz] = curveGen(Fpz, Oz, T7, plane, surf, dt);
    %display(['Fpz-T7-Oz curve length: ' num2str(len_FpzT7Oz)]);

    [Fp1h,~] = curveWalk(curve_pts_FpzT7Oz, 0.5*len_FpzT7Oz/10);
    Fp1h = bringPtsToSurf(surf,Fp1h);
    %display(['Fp1h = [' num2str(Fp1h(1,1)) ' ' num2str(Fp1h(1,2)) ' ' num2str(Fp1h(1,3)) '], ' num2str(len) ' away from Fpz']);
    
    [Fp1,~] = curveWalk(curve_pts_FpzT7Oz, 1*len_FpzT7Oz/10);
    Fp1 = bringPtsToSurf(surf,Fp1);
    %display(['Fp1 = [' num2str(Fp1(1,1)) ' ' num2str(Fp1(1,2)) ' ' num2str(Fp1(1,3)) '], ' num2str(len) ' away from Fpz']);

    [AFp7,~] = curveWalk(curve_pts_FpzT7Oz, 1.5*len_FpzT7Oz/10);
    AFp7 = bringPtsToSurf(surf,AFp7);
    %display(['AFp7 = [' num2str(AFp7(1,1)) ' ' num2str(AFp7(1,2)) ' ' num2str(AFp7(1,3)) '], ' num2str(len) ' away from Fpz']);
    
    [AF7,~] = curveWalk(curve_pts_FpzT7Oz, 2*len_FpzT7Oz/10);
    AF7 = bringPtsToSurf(surf,AF7);
    %display(['AF7 = [' num2str(AF7(1,1)) ' ' num2str(AF7(1,2)) ' ' num2str(AF7(1,3)) '], ' num2str(len) ' away from Fpz']);
    
    [AFF7,~] = curveWalk(curve_pts_FpzT7Oz, 2.5*len_FpzT7Oz/10);
    AFF7 = bringPtsToSurf(surf,AFF7);
    %display(['AFF7 = [' num2str(AFF7(1,1)) ' ' num2str(AFF7(1,2)) ' ' num2str(AFF7(1,3)) '], ' num2str(len) ' away from Fpz']);
    
    [F7,~] = curveWalk(curve_pts_FpzT7Oz, 3*len_FpzT7Oz/10);
    F7 = bringPtsToSurf(surf,F7);
    %display(['F7 = [' num2str(F7(1,1)) ' ' num2str(F7(1,2)) ' ' num2str(F7(1,3)) '], ' num2str(len) ' away from Fpz']);
    
    [FFT7,~] = curveWalk(curve_pts_FpzT7Oz, 3.5*len_FpzT7Oz/10);
    FFT7 = bringPtsToSurf(surf,FFT7);
    %display(['FFT7 = [' num2str(FFT7(1,1)) ' ' num2str(FFT7(1,2)) ' ' num2str(FFT7(1,3)) '], ' num2str(len) ' away from Fpz']);

    [FT7,~] = curveWalk(curve_pts_FpzT7Oz, 4*len_FpzT7Oz/10);
    FT7 = bringPtsToSurf(surf,FT7);
    %display(['FT7 = [' num2str(FT7(1,1)) ' ' num2str(FT7(1,2)) ' ' num2str(FT7(1,3)) '], ' num2str(len) ' away from Fpz']);
    
    [FTT7,~] = curveWalk(curve_pts_FpzT7Oz, 4.5*len_FpzT7Oz/10);
    FTT7 = bringPtsToSurf(surf,FTT7);
    %display(['FTT7 = [' num2str(FTT7(1,1)) ' ' num2str(FTT7(1,2)) ' ' num2str(FTT7(1,3)) '], ' num2str(len) ' away from Fpz']);
    
    [TTP7,~] = curveWalk(curve_pts_FpzT7Oz, 5.5*len_FpzT7Oz/10);
    TTP7 = bringPtsToSurf(surf,TTP7);
    %display(['TTP7 = [' num2str(TTP7(1,1)) ' ' num2str(TTP7(1,2)) ' ' num2str(TTP7(1,3)) '], ' num2str(len) ' away from Fpz']);
    
    [TP7,~] = curveWalk(curve_pts_FpzT7Oz, 6*len_FpzT7Oz/10);
    TP7 = bringPtsToSurf(surf,TP7);
    %display(['TP7 = [' num2str(TP7(1,1)) ' ' num2str(TP7(1,2)) ' ' num2str(TP7(1,3)) '], ' num2str(len) ' away from Fpz']);
    
    [TPP7,~] = curveWalk(curve_pts_FpzT7Oz, 6.5*len_FpzT7Oz/10);
    TPP7 = bringPtsToSurf(surf,TPP7);
    %display(['TPP7 = [' num2str(TPP7(1,1)) ' ' num2str(TPP7(1,2)) ' ' num2str(TPP7(1,3)) '], ' num2str(len) ' away from Fpz']);
    
    [P7,~] = curveWalk(curve_pts_FpzT7Oz, 7*len_FpzT7Oz/10);
    P7 = bringPtsToSurf(surf,P7);
    %display(['P7 = [' num2str(P7(1,1)) ' ' num2str(P7(1,2)) ' ' num2str(P7(1,3)) '], ' num2str(len) ' away from Fpz']);
    
    [PPO7,~] = curveWalk(curve_pts_FpzT7Oz, 7.5*len_FpzT7Oz/10);
    PPO7 = bringPtsToSurf(surf,PPO7);
    %display(['PPO7 = [' num2str(PPO7(1,1)) ' ' num2str(PPO7(1,2)) ' ' num2str(PPO7(1,3)) '], ' num2str(len) ' away from Fpz']);
    
    [PO7,~] = curveWalk(curve_pts_FpzT7Oz, 8*len_FpzT7Oz/10);
    PO7 = bringPtsToSurf(surf,PO7);
    %display(['PO7 = [' num2str(PO7(1,1)) ' ' num2str(PO7(1,2)) ' ' num2str(PO7(1,3)) '], ' num2str(len) ' away from Fpz']);
    
    [POO7,~] = curveWalk(curve_pts_FpzT7Oz, 8.5*len_FpzT7Oz/10);
    POO7 = bringPtsToSurf(surf,POO7);
    %display(['POO7 = [' num2str(POO7(1,1)) ' ' num2str(POO7(1,2)) ' ' num2str(POO7(1,3)) '], ' num2str(len) ' away from Fpz']);

    [O1,~] = curveWalk(curve_pts_FpzT7Oz, 9*len_FpzT7Oz/10);
    O1 = bringPtsToSurf(surf,O1);
    %display(['O1 = [' num2str(O1(1,1)) ' ' num2str(O1(1,2)) ' ' num2str(O1(1,3)) '], ' num2str(len) ' away from Fpz']);
    
    [O1h,~] = curveWalk(curve_pts_FpzT7Oz, 9.5*len_FpzT7Oz/10);
    O1h = bringPtsToSurf(surf,O1h);
    %display(['O1h = [' num2str(O1h(1,1)) ' ' num2str(O1h(1,2)) ' ' num2str(O1h(1,3)) '], ' num2str(len) ' away from Fpz']);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find point along Fpz-T8-Oz curve 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    plane = plane_equation(Fpz, T8, Oz);
    [curve_pts_FpzT8Oz,len_FpzT8Oz] = curveGen(Fpz, Oz, T8, plane, surf, dt);
    %display(['Fpz-T8-Oz curve length: ' num2str(len_FpzT8Oz)]);

    [Fp2h,~] = curveWalk(curve_pts_FpzT8Oz, 0.5*len_FpzT8Oz/10);
    Fp2h = bringPtsToSurf(surf,Fp2h);
    %display(['Fp2h = [' num2str(Fp2h(1,1)) ' ' num2str(Fp2h(1,2)) ' ' num2str(Fp2h(1,3)) '], ' num2str(len) ' away from Fpz']);
    
    [Fp2,~] = curveWalk(curve_pts_FpzT8Oz, 1*len_FpzT8Oz/10);
    Fp2 = bringPtsToSurf(surf,Fp2);
    %display(['Fp2 = [' num2str(Fp2(1,1)) ' ' num2str(Fp2(1,2)) ' ' num2str(Fp2(1,3)) '], ' num2str(len) ' away from Fpz']);

    [AFp8,~] = curveWalk(curve_pts_FpzT8Oz, 1.5*len_FpzT8Oz/10);
    AFp8 = bringPtsToSurf(surf,AFp8);
    %display(['AFp8 = [' num2str(AFp8(1,1)) ' ' num2str(AFp8(1,2)) ' ' num2str(AFp8(1,3)) '], ' num2str(len) ' away from Fpz']);
    
    [AF8,~] = curveWalk(curve_pts_FpzT8Oz, 2*len_FpzT8Oz/10);
    AF8 = bringPtsToSurf(surf,AF8);
    %display(['AF8 = [' num2str(AF8(1,1)) ' ' num2str(AF8(1,2)) ' ' num2str(AF8(1,3)) '], ' num2str(len) ' away from Fpz']);
    
    [AFF8,~] = curveWalk(curve_pts_FpzT8Oz, 2.5*len_FpzT8Oz/10);
    AFF8 = bringPtsToSurf(surf,AFF8);
    %display(['AFF8 = [' num2str(AFF8(1,1)) ' ' num2str(AFF8(1,2)) ' ' num2str(AFF8(1,3)) '], ' num2str(len) ' away from Fpz']);
    
    [F8,~] = curveWalk(curve_pts_FpzT8Oz, 3*len_FpzT8Oz/10);
    F8 = bringPtsToSurf(surf,F8);
    %display(['F8 = [' num2str(F8(1,1)) ' ' num2str(F8(1,2)) ' ' num2str(F8(1,3)) '], ' num2str(len) ' away from Fpz']);
    
    [FFT8,~] = curveWalk(curve_pts_FpzT8Oz, 3.5*len_FpzT8Oz/10);
    FFT8 = bringPtsToSurf(surf,FFT8);
    %display(['FFT8 = [' num2str(FFT8(1,1)) ' ' num2str(FFT8(1,2)) ' ' num2str(FFT8(1,3)) '], ' num2str(len) ' away from Fpz']);

    [FT8,~] = curveWalk(curve_pts_FpzT8Oz, 4*len_FpzT8Oz/10);
    FT8 = bringPtsToSurf(surf,FT8);
    %display(['FT8 = [' num2str(FT8(1,1)) ' ' num2str(FT8(1,2)) ' ' num2str(FT8(1,3)) '], ' num2str(len) ' away from Fpz']);
    
    [FTT8,~] = curveWalk(curve_pts_FpzT8Oz, 4.5*len_FpzT8Oz/10);
    FTT8 = bringPtsToSurf(surf,FTT8);
    %display(['FTT8 = [' num2str(FTT8(1,1)) ' ' num2str(FTT8(1,2)) ' ' num2str(FTT8(1,3)) '], ' num2str(len) ' away from Fpz']);
    
    [TTP8,~] = curveWalk(curve_pts_FpzT8Oz, 5.5*len_FpzT8Oz/10);
    TTP8 = bringPtsToSurf(surf,TTP8);
    %display(['TTP8 = [' num2str(TTP8(1,1)) ' ' num2str(TTP8(1,2)) ' ' num2str(TTP8(1,3)) '], ' num2str(len) ' away from Fpz']);
    
    [TP8,~] = curveWalk(curve_pts_FpzT8Oz, 6*len_FpzT8Oz/10);
    TP8 = bringPtsToSurf(surf,TP8);
    %display(['TP8 = [' num2str(TP8(1,1)) ' ' num2str(TP8(1,2)) ' ' num2str(TP8(1,3)) '], ' num2str(len) ' away from Fpz']);
    
    [TPP8,~] = curveWalk(curve_pts_FpzT8Oz, 6.5*len_FpzT8Oz/10);
    TPP8 = bringPtsToSurf(surf,TPP8);
    %display(['TPP8 = [' num2str(TPP8(1,1)) ' ' num2str(TPP8(1,2)) ' ' num2str(TPP8(1,3)) '], ' num2str(len) ' away from Fpz']);
    
    [P8,~] = curveWalk(curve_pts_FpzT8Oz, 7*len_FpzT8Oz/10);
    P8 = bringPtsToSurf(surf,P8);
    %display(['P8 = [' num2str(P8(1,1)) ' ' num2str(P8(1,2)) ' ' num2str(P8(1,3)) '], ' num2str(len) ' away from Fpz']);
    
    [PPO8,~] = curveWalk(curve_pts_FpzT8Oz, 7.5*len_FpzT8Oz/10);
    PPO8 = bringPtsToSurf(surf,PPO8);
    %display(['PPO8 = [' num2str(PPO8(1,1)) ' ' num2str(PPO8(1,2)) ' ' num2str(PPO8(1,3)) '], ' num2str(len) ' away from Fpz']);
    
    [PO8,~] = curveWalk(curve_pts_FpzT8Oz, 8*len_FpzT8Oz/10);
    PO8 = bringPtsToSurf(surf,PO8);
    %display(['PO8 = [' num2str(PO8(1,1)) ' ' num2str(PO8(1,2)) ' ' num2str(PO8(1,3)) '], ' num2str(len) ' away from Fpz']);
    
    [POO8,~] = curveWalk(curve_pts_FpzT8Oz, 8.5*len_FpzT8Oz/10);
    POO8 = bringPtsToSurf(surf,POO8);
    %display(['POO8 = [' num2str(POO8(1,1)) ' ' num2str(POO8(1,2)) ' ' num2str(POO8(1,3)) '], ' num2str(len) ' away from Fpz']);

    [O2,~] = curveWalk(curve_pts_FpzT8Oz, 9*len_FpzT8Oz/10);
    O2 = bringPtsToSurf(surf,O2);
    %display(['O2 = [' num2str(O2(1,1)) ' ' num2str(O2(1,2)) ' ' num2str(O2(1,3)) '], ' num2str(len) ' away from Fpz']);
    
    [O2h,~] = curveWalk(curve_pts_FpzT8Oz, 9.5*len_FpzT8Oz/10);
    O2h = bringPtsToSurf(surf,O2h);
    %display(['O2h = [' num2str(O2h(1,1)) ' ' num2str(O2h(1,2)) ' ' num2str(O2h(1,3)) '], ' num2str(len) ' away from Fpz']);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find point along F7-Fz-F8 curve 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    plane = plane_equation(F7, Fz, F8);
    [curve_pts_F7FzF8,len_F7FzF8] = curveGen(F7, F8, Fz, plane, surf);
    %display(['F7-Fz-F8 curve length: ' num2str(len_F7FzF8)]);
    
    [F7h,~] = curveWalk(curve_pts_F7FzF8, 0.5*len_F7FzF8/8);
    F7h = bringPtsToSurf(surf,F7h);
    %display(['F7h = [' num2str(F7h(1,1)) ' ' num2str(F7h(1,2)) ' ' num2str(F7h(1,3)) '], ' num2str(len) ' away from F7']);

    [F5,~] = curveWalk(curve_pts_F7FzF8, 1*len_F7FzF8/8);
    F5 = bringPtsToSurf(surf,F5);
    %display(['F5 = [' num2str(F5(1,1)) ' ' num2str(F5(1,2)) ' ' num2str(F5(1,3)) '], ' num2str(len) ' away from F7']);
    
    [F5h,~] = curveWalk(curve_pts_F7FzF8, 1.5*len_F7FzF8/8);
    F5h = bringPtsToSurf(surf,F5h);
    %display(['F5h = [' num2str(F5h(1,1)) ' ' num2str(F5h(1,2)) ' ' num2str(F5h(1,3)) '], ' num2str(len) ' away from F7']);
    
    [F3,~] = curveWalk(curve_pts_F7FzF8, 2*len_F7FzF8/8);
    F3 = bringPtsToSurf(surf,F3);
    %display(['F3 = [' num2str(F3(1,1)) ' ' num2str(F3(1,2)) ' ' num2str(F3(1,3)) '], ' num2str(len) ' away from F7']);

    [F3h,~] = curveWalk(curve_pts_F7FzF8, 2.5*len_F7FzF8/8);
    F3h = bringPtsToSurf(surf,F3h);
    %display(['F3h = [' num2str(F3h(1,1)) ' ' num2str(F3h(1,2)) ' ' num2str(F3h(1,3)) '], ' num2str(len) ' away from F7']);
    
    [F1,~] = curveWalk(curve_pts_F7FzF8, 3*len_F7FzF8/8);
    F1 = bringPtsToSurf(surf,F1);
    %display(['F1 = [' num2str(F1(1,1)) ' ' num2str(F1(1,2)) ' ' num2str(F1(1,3)) '], ' num2str(len) ' away from F7']);
    
    [F1h,~] = curveWalk(curve_pts_F7FzF8, 3.5*len_F7FzF8/8);
    F1h = bringPtsToSurf(surf,F1h);
    %display(['F1h = [' num2str(F1h(1,1)) ' ' num2str(F1h(1,2)) ' ' num2str(F1h(1,3)) '], ' num2str(len) ' away from F7']);
    
    [F2h,~] = curveWalk(curve_pts_F7FzF8, 4.5*len_F7FzF8/8);
    F2h = bringPtsToSurf(surf,F2h);
    %display(['F2h = [' num2str(F2h(1,1)) ' ' num2str(F2h(1,2)) ' ' num2str(F2h(1,3)) '], ' num2str(len) ' away from F8']);
    
    [F2,~] = curveWalk(curve_pts_F7FzF8, 5*len_F7FzF8/8);
    F2 = bringPtsToSurf(surf,F2);
    %display(['F2 = [' num2str(F2(1,1)) ' ' num2str(F2(1,2)) ' ' num2str(F2(1,3)) '], ' num2str(len) ' away from F8']);
    
    [F4h,~] = curveWalk(curve_pts_F7FzF8, 5.5*len_F7FzF8/8);
    F4h = bringPtsToSurf(surf,F4h);
    %display(['F4h = [' num2str(F4h(1,1)) ' ' num2str(F4h(1,2)) ' ' num2str(F4h(1,3)) '], ' num2str(len) ' away from F8']);
    
    [F4,~] = curveWalk(curve_pts_F7FzF8, 6*len_F7FzF8/8);
    F4 = bringPtsToSurf(surf,F4);
    %display(['F4 = [' num2str(F4(1,1)) ' ' num2str(F4(1,2)) ' ' num2str(F4(1,3)) '], ' num2str(len) ' away from F8']);
    
    [F6h,~] = curveWalk(curve_pts_F7FzF8, 6.5*len_F7FzF8/8);
    F6h = bringPtsToSurf(surf,F6h);
    %display(['F6h = [' num2str(F6h(1,1)) ' ' num2str(F6h(1,2)) ' ' num2str(F6h(1,3)) '], ' num2str(len) ' away from F8']);
    
    [F6,~] = curveWalk(curve_pts_F7FzF8, 7*len_F7FzF8/8);
    F6 = bringPtsToSurf(surf,F6);
    %display(['F6 = [' num2str(F6(1,1)) ' ' num2str(F6(1,2)) ' ' num2str(F6(1,3)) '], ' num2str(len) ' away from F8']);
    
    [F8h,~] = curveWalk(curve_pts_F7FzF8, 7.5*len_F7FzF8/8);
    F8h = bringPtsToSurf(surf,F8h);
    %display(['F8h = [' num2str(F8h(1,1)) ' ' num2str(F8h(1,2)) ' ' num2str(F8h(1,3)) '], ' num2str(len) ' away from F8']);
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find point along P7-Pz-P8 curve 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    plane = plane_equation(P7, Pz, P8);
    [curve_pts_P7PzP8,len_P7PzP8] = curveGen(P7, P8, Pz, plane, surf,dt);
    %display(['P7-Pz-P8 curve length: ' num2str(len_P7PzP8)]);

    [P7h,~] = curveWalk(curve_pts_P7PzP8, 0.5*len_P7PzP8/8);
    P7h = bringPtsToSurf(surf,P7h);
    %display(['P7h = [' num2str(P7h(1,1)) ' ' num2str(P7h(1,2)) ' ' num2str(P7h(1,3)) '], ' num2str(len) ' away from P7']);
    
    [P5,~] = curveWalk(curve_pts_P7PzP8, 1*len_P7PzP8/8);
    P5 = bringPtsToSurf(surf,P5);
    %display(['P5 = [' num2str(P5(1,1)) ' ' num2str(P5(1,2)) ' ' num2str(P5(1,3)) '], ' num2str(len) ' away from P7']);
    
    [P5h,~] = curveWalk(curve_pts_P7PzP8, 1.5*len_P7PzP8/8);
    P5h = bringPtsToSurf(surf,P5h);
    %display(['P5h = [' num2str(P5h(1,1)) ' ' num2str(P5h(1,2)) ' ' num2str(P5h(1,3)) '], ' num2str(len) ' away from P7']);
    
    [P3,~] = curveWalk(curve_pts_P7PzP8, 2*len_P7PzP8/8);
    P3 = bringPtsToSurf(surf,P3);
    %display(['P3 = [' num2str(P3(1,1)) ' ' num2str(P3(1,2)) ' ' num2str(P3(1,3)) '], ' num2str(len) ' away from P7']);

    [P3h,~] = curveWalk(curve_pts_P7PzP8, 2.5*len_P7PzP8/8);
    P3h = bringPtsToSurf(surf,P3h);
    %display(['P3h = [' num2str(P3h(1,1)) ' ' num2str(P3h(1,2)) ' ' num2str(P3h(1,3)) '], ' num2str(len) ' away from P7']);
    
    [P1,~] = curveWalk(curve_pts_P7PzP8, 3*len_P7PzP8/8);
    P1 = bringPtsToSurf(surf,P1);
    %display(['P1 = [' num2str(P1(1,1)) ' ' num2str(P1(1,2)) ' ' num2str(P1(1,3)) '], ' num2str(len) ' away from P7']);
    
    [P1h,~] = curveWalk(curve_pts_P7PzP8, 3.5*len_P7PzP8/8);
    P1h = bringPtsToSurf(surf,P1h);
    %display(['P1h = [' num2str(P1h(1,1)) ' ' num2str(P1h(1,2)) ' ' num2str(P1h(1,3)) '], ' num2str(len) ' away from P7']);
    
    [P2h,~] = curveWalk(curve_pts_P7PzP8, 4.5*len_P7PzP8/8);
    P2h = bringPtsToSurf(surf,P2h);
    %display(['P2h = [' num2str(P2h(1,1)) ' ' num2str(P2h(1,2)) ' ' num2str(P2h(1,3)) '], ' num2str(len) ' away from P8']);
    
    [P2,~] = curveWalk(curve_pts_P7PzP8, 5*len_P7PzP8/8);
    P2 = bringPtsToSurf(surf,P2);
    %display(['P2 = [' num2str(P2(1,1)) ' ' num2str(P2(1,2)) ' ' num2str(P2(1,3)) '], ' num2str(len) ' away from P8']);
    
    [P4h,~] = curveWalk(curve_pts_P7PzP8, 5.5*len_P7PzP8/8);
    P4h = bringPtsToSurf(surf,P4h);
    %display(['P4h = [' num2str(P4h(1,1)) ' ' num2str(P4h(1,2)) ' ' num2str(P4h(1,3)) '], ' num2str(len) ' away from P8']);
    
    [P4,~] = curveWalk(curve_pts_P7PzP8, 6*len_P7PzP8/8);
    P4 = bringPtsToSurf(surf,P4);
    %display(['P4 = [' num2str(P4(1,1)) ' ' num2str(P4(1,2)) ' ' num2str(P4(1,3)) '], ' num2str(len) ' away from P8']);
    
    [P6h,~] = curveWalk(curve_pts_P7PzP8, 6.5*len_P7PzP8/8);
    P6h = bringPtsToSurf(surf,P6h);
    %display(['P6h = [' num2str(P6h(1,1)) ' ' num2str(P6h(1,2)) ' ' num2str(P6h(1,3)) '], ' num2str(len) ' away from P8']);
    
    [P6,~] = curveWalk(curve_pts_P7PzP8, 7*len_P7PzP8/8);
    P6 = bringPtsToSurf(surf,P6);
    %display(['P6 = [' num2str(P6(1,1)) ' ' num2str(P6(1,2)) ' ' num2str(P6(1,3)) '], ' num2str(len) ' away from P8']);
    
    [P8h,~] = curveWalk(curve_pts_P7PzP8, 7.5*len_P7PzP8/8);
    P8h = bringPtsToSurf(surf,P8h);
    %display(['P8h = [' num2str(P8h(1,1)) ' ' num2str(P8h(1,2)) ' ' num2str(P8h(1,3)) '], ' num2str(len) ' away from P8']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find point along AF7-PFz-AF8 curve 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    plane = plane_equation(AF7, AFz, AF8);
    [curve_pts_AF7AFzAF8,len_AF7AFzAF8] = curveGen(AF7, AF8, AFz, plane, surf,dt);
    %display(['AF7-AFz-AF8 curve length: ' num2str(len_AF7AFzAF8)]);
    
    [AF7h,~] = curveWalk(curve_pts_AF7AFzAF8, 0.5*len_AF7AFzAF8/8);
    AF7h = bringPtsToSurf(surf,AF7h);
    %display(['AF7h = [' num2str(AF7h(1,1)) ' ' num2str(AF7h(1,2)) ' ' num2str(AF7h(1,3)) '], ' num2str(len) ' away from AF7']);

    [AF5,~] = curveWalk(curve_pts_AF7AFzAF8, 1*len_AF7AFzAF8/8);
    AF5 = bringPtsToSurf(surf,AF5);
    %display(['AF5 = [' num2str(AF5(1,1)) ' ' num2str(AF5(1,2)) ' ' num2str(AF5(1,3)) '], ' num2str(len) ' away from AF7']);
    
    [AF5h,~] = curveWalk(curve_pts_AF7AFzAF8, 1.5*len_AF7AFzAF8/8);
    AF5h = bringPtsToSurf(surf,AF5h);
    %display(['AF5h = [' num2str(AF5h(1,1)) ' ' num2str(AF5h(1,2)) ' ' num2str(AF5h(1,3)) '], ' num2str(len) ' away from AF7']);
    
    [AF3,~] = curveWalk(curve_pts_AF7AFzAF8, 2*len_AF7AFzAF8/8);
    AF3 = bringPtsToSurf(surf,AF3);
    %display(['AF3 = [' num2str(AF3(1,1)) ' ' num2str(AF3(1,2)) ' ' num2str(AF3(1,3)) '], ' num2str(len) ' away from AF7']);

    [AF3h,~] = curveWalk(curve_pts_AF7AFzAF8, 2.5*len_AF7AFzAF8/8);
    AF3h = bringPtsToSurf(surf,AF3h);
    %display(['AF3h = [' num2str(AF3h(1,1)) ' ' num2str(AF3h(1,2)) ' ' num2str(AF3h(1,3)) '], ' num2str(len) ' away from AF7']);
    
    [AF1,~] = curveWalk(curve_pts_AF7AFzAF8, 3*len_AF7AFzAF8/8);
    AF1 = bringPtsToSurf(surf,AF1);
    %display(['AF1 = [' num2str(AF1(1,1)) ' ' num2str(AF1(1,2)) ' ' num2str(AF1(1,3)) '], ' num2str(len) ' away from AF7']);
    
    [AF1h,~] = curveWalk(curve_pts_AF7AFzAF8, 3.5*len_AF7AFzAF8/8);
    AF1h = bringPtsToSurf(surf,AF1h);
    %display(['AF1h = [' num2str(AF1h(1,1)) ' ' num2str(AF1h(1,2)) ' ' num2str(AF1h(1,3)) '], ' num2str(len) ' away from AF7']);
    
    [AF2h,~] = curveWalk(curve_pts_AF7AFzAF8, 4.5*len_AF7AFzAF8/8);
    AF2h = bringPtsToSurf(surf,AF2h);
    %display(['AF2h = [' num2str(AF2h(1,1)) ' ' num2str(AF2h(1,2)) ' ' num2str(AF2h(1,3)) '], ' num2str(len) ' away from AF8']);
    
    [AF2,~] = curveWalk(curve_pts_AF7AFzAF8, 5*len_AF7AFzAF8/8);
    AF2 = bringPtsToSurf(surf,AF2);
    %display(['AF2 = [' num2str(AF2(1,1)) ' ' num2str(AF2(1,2)) ' ' num2str(AF2(1,3)) '], ' num2str(len) ' away from AF8']);
    
    [AF4h,~] = curveWalk(curve_pts_AF7AFzAF8, 5.5*len_AF7AFzAF8/8);
    AF4h = bringPtsToSurf(surf,AF4h);
    %display(['AF4h = [' num2str(AF4h(1,1)) ' ' num2str(AF4h(1,2)) ' ' num2str(AF4h(1,3)) '], ' num2str(len) ' away from AF8']);
    
    [AF4,~] = curveWalk(curve_pts_AF7AFzAF8, 6*len_AF7AFzAF8/8);
    AF4 = bringPtsToSurf(surf,AF4);
    %display(['AF4 = [' num2str(AF4(1,1)) ' ' num2str(AF4(1,2)) ' ' num2str(AF4(1,3)) '], ' num2str(len) ' away from AF8']);
    
    [AF6h,~] = curveWalk(curve_pts_AF7AFzAF8, 6.5*len_AF7AFzAF8/8);
    AF6h = bringPtsToSurf(surf,AF6h);
    %display(['AF6h = [' num2str(AF6h(1,1)) ' ' num2str(AF6h(1,2)) ' ' num2str(AF6h(1,3)) '], ' num2str(len) ' away from AF8']);
    
    [AF6,~] = curveWalk(curve_pts_AF7AFzAF8, 7*len_AF7AFzAF8/8);
    AF6 = bringPtsToSurf(surf,AF6);
    %display(['AF6 = [' num2str(AF6(1,1)) ' ' num2str(AF6(1,2)) ' ' num2str(AF6(1,3)) '], ' num2str(len) ' away from AF8']);
    
    [AF8h,~] = curveWalk(curve_pts_AF7AFzAF8, 7.5*len_AF7AFzAF8/8);
    AF8h = bringPtsToSurf(surf,AF8h);
    %display(['AF8h = [' num2str(AF8h(1,1)) ' ' num2str(AF8h(1,2)) ' ' num2str(AF8h(1,3)) '], ' num2str(len) ' away from AF8']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find point along PO7-POz-PO8 curve 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    plane = plane_equation(PO7, POz, PO8);
    [curve_pts_PO7POzPO8,len_PO7POzPO8] = curveGen(PO7, PO8, POz, plane, surf,dt);
    %display(['PO7-POz-PO8 curve length: ' num2str(len_PO7POzPO8)]);

    [PO7h,~] = curveWalk(curve_pts_PO7POzPO8, 0.5*len_PO7POzPO8/8);
    PO7h = bringPtsToSurf(surf,PO7h);
    %display(['PO7h = [' num2str(PO7h(1,1)) ' ' num2str(PO7h(1,2)) ' ' num2str(PO7h(1,3)) '], ' num2str(len) ' away from PO7']);
    
    [PO5,~] = curveWalk(curve_pts_PO7POzPO8, 1*len_PO7POzPO8/8);
    PO5 = bringPtsToSurf(surf,PO5);
    %display(['PO5 = [' num2str(PO5(1,1)) ' ' num2str(PO5(1,2)) ' ' num2str(PO5(1,3)) '], ' num2str(len) ' away from PO7']);
    
    [PO5h,~] = curveWalk(curve_pts_PO7POzPO8, 1.5*len_PO7POzPO8/8);
    PO5h = bringPtsToSurf(surf,PO5h);
    %display(['PO5h = [' num2str(PO5h(1,1)) ' ' num2str(PO5h(1,2)) ' ' num2str(PO5h(1,3)) '], ' num2str(len) ' away from PO7']);
    
    [PO3,~] = curveWalk(curve_pts_PO7POzPO8, 2*len_PO7POzPO8/8);
    PO3 = bringPtsToSurf(surf,PO3);
    %display(['PO3 = [' num2str(PO3(1,1)) ' ' num2str(PO3(1,2)) ' ' num2str(PO3(1,3)) '], ' num2str(len) ' away from PO7']);

    [PO3h,~] = curveWalk(curve_pts_PO7POzPO8, 2.5*len_PO7POzPO8/8);
    PO3h = bringPtsToSurf(surf,PO3h);
    %display(['PO3h = [' num2str(PO3h(1,1)) ' ' num2str(PO3h(1,2)) ' ' num2str(PO3h(1,3)) '], ' num2str(len) ' away from PO7']);
    
    [PO1,~] = curveWalk(curve_pts_PO7POzPO8, 3*len_PO7POzPO8/8);
    PO1 = bringPtsToSurf(surf,PO1);
    %display(['PO1 = [' num2str(PO1(1,1)) ' ' num2str(PO1(1,2)) ' ' num2str(PO1(1,3)) '], ' num2str(len) ' away from PO7']);
    
    [PO1h,~] = curveWalk(curve_pts_PO7POzPO8, 3.5*len_PO7POzPO8/8);
    PO1h = bringPtsToSurf(surf,PO1h);
    %display(['PO1h = [' num2str(PO1h(1,1)) ' ' num2str(PO1h(1,2)) ' ' num2str(PO1h(1,3)) '], ' num2str(len) ' away from PO7']);
    
    [PO2h,~] = curveWalk(curve_pts_PO7POzPO8, 4.5*len_PO7POzPO8/8);
    PO2h = bringPtsToSurf(surf,PO2h);
    %display(['PO2h = [' num2str(PO2h(1,1)) ' ' num2str(PO2h(1,2)) ' ' num2str(PO2h(1,3)) '], ' num2str(len) ' away from PO8']);
    
    [PO2,~] = curveWalk(curve_pts_PO7POzPO8, 5*len_PO7POzPO8/8);
    PO2 = bringPtsToSurf(surf,PO2);
    %display(['PO2 = [' num2str(PO2(1,1)) ' ' num2str(PO2(1,2)) ' ' num2str(PO2(1,3)) '], ' num2str(len) ' away from PO8']);
    
    [PO4h,~] = curveWalk(curve_pts_PO7POzPO8, 5.5*len_PO7POzPO8/8);
    PO4h = bringPtsToSurf(surf,PO4h);
    %display(['PO4h = [' num2str(PO4h(1,1)) ' ' num2str(PO4h(1,2)) ' ' num2str(PO4h(1,3)) '], ' num2str(len) ' away from PO8']);
    
    [PO4,~] = curveWalk(curve_pts_PO7POzPO8, 6*len_PO7POzPO8/8);
    PO4 = bringPtsToSurf(surf,PO4);
    %display(['PO4 = [' num2str(PO4(1,1)) ' ' num2str(PO4(1,2)) ' ' num2str(PO4(1,3)) '], ' num2str(len) ' away from PO8']);
    
    [PO6h,~] = curveWalk(curve_pts_PO7POzPO8, 6.5*len_PO7POzPO8/8);
    PO6h = bringPtsToSurf(surf,PO6h);
    %display(['PO6h = [' num2str(PO6h(1,1)) ' ' num2str(PO6h(1,2)) ' ' num2str(PO6h(1,3)) '], ' num2str(len) ' away from PO8']);
    
    [PO6,~] = curveWalk(curve_pts_PO7POzPO8, 7*len_PO7POzPO8/8);
    PO6 = bringPtsToSurf(surf,PO6);
    %display(['PO6 = [' num2str(PO6(1,1)) ' ' num2str(PO6(1,2)) ' ' num2str(PO6(1,3)) '], ' num2str(len) ' away from PO8']);
    
    [PO8h,~] = curveWalk(curve_pts_PO7POzPO8, 7.5*len_PO7POzPO8/8);
    PO8h = bringPtsToSurf(surf,PO8h);
    %display(['PO8h = [' num2str(PO8h(1,1)) ' ' num2str(PO8h(1,2)) ' ' num2str(PO8h(1,3)) '], ' num2str(len) ' away from PO8']);
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find point along FT7-FCz-FT8 curve 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    plane = plane_equation(FT7, FCz, FT8);
    [curve_pts_FT7FCzFT8,len_FT7FCzFT8] = curveGen(FT7, FT8, FCz, plane, surf,dt);
    %display(['FT7-FCz-FT8 curve length: ' num2str(len_FT7FCzFT8)]);
    
    [FT7h,~] = curveWalk(curve_pts_FT7FCzFT8, 0.5*len_FT7FCzFT8/8);
    FT7h = bringPtsToSurf(surf,FT7h);
    %display(['FT7h = [' num2str(FT7h(1,1)) ' ' num2str(FT7h(1,2)) ' ' num2str(FT7h(1,3)) '], ' num2str(len) ' away from FT7']);

    [FC5,~] = curveWalk(curve_pts_FT7FCzFT8, 1*len_FT7FCzFT8/8);
    FC5 = bringPtsToSurf(surf,FC5);
    %display(['FC5 = [' num2str(FC5(1,1)) ' ' num2str(FC5(1,2)) ' ' num2str(FC5(1,3)) '], ' num2str(len) ' away from FT7']);
    
    [FC5h,~] = curveWalk(curve_pts_FT7FCzFT8, 1.5*len_FT7FCzFT8/8);
    FC5h = bringPtsToSurf(surf,FC5h);
    %display(['FC5h = [' num2str(FC5h(1,1)) ' ' num2str(FC5h(1,2)) ' ' num2str(FC5h(1,3)) '], ' num2str(len) ' away from FT7']);
    
    [FC3,~] = curveWalk(curve_pts_FT7FCzFT8, 2*len_FT7FCzFT8/8);
    FC3 = bringPtsToSurf(surf,FC3);
    %display(['FC3 = [' num2str(FC3(1,1)) ' ' num2str(FC3(1,2)) ' ' num2str(FC3(1,3)) '], ' num2str(len) ' away from FT7']);

    [FC3h,~] = curveWalk(curve_pts_FT7FCzFT8, 2.5*len_FT7FCzFT8/8);
    FC3h = bringPtsToSurf(surf,FC3h);
    %display(['FC3h = [' num2str(FC3h(1,1)) ' ' num2str(FC3h(1,2)) ' ' num2str(FC3h(1,3)) '], ' num2str(len) ' away from FT7']);
    
    [FC1,~] = curveWalk(curve_pts_FT7FCzFT8, 3*len_FT7FCzFT8/8);
    FC1 = bringPtsToSurf(surf,FC1);
    %display(['FC1 = [' num2str(FC1(1,1)) ' ' num2str(FC1(1,2)) ' ' num2str(FC1(1,3)) '], ' num2str(len) ' away from FT7']);
    
    [FC1h,~] = curveWalk(curve_pts_FT7FCzFT8, 3.5*len_FT7FCzFT8/8);
    FC1h = bringPtsToSurf(surf,FC1h);
    %display(['FC1h = [' num2str(FC1h(1,1)) ' ' num2str(FC1h(1,2)) ' ' num2str(FC1h(1,3)) '], ' num2str(len) ' away from FT7']);
    
    [FC2h,~] = curveWalk(curve_pts_FT7FCzFT8, 4.5*len_FT7FCzFT8/8);
    FC2h = bringPtsToSurf(surf,FC2h);
    %display(['FC2h = [' num2str(FC2h(1,1)) ' ' num2str(FC2h(1,2)) ' ' num2str(FC2h(1,3)) '], ' num2str(len) ' away from FT8']);
    
    [FC2,~] = curveWalk(curve_pts_FT7FCzFT8, 5*len_FT7FCzFT8/8);
    FC2 = bringPtsToSurf(surf,FC2);
    %display(['FC2 = [' num2str(FC2(1,1)) ' ' num2str(FC2(1,2)) ' ' num2str(FC2(1,3)) '], ' num2str(len) ' away from FT8']);
    
    [FC4h,~] = curveWalk(curve_pts_FT7FCzFT8, 5.5*len_FT7FCzFT8/8);
    FC4h = bringPtsToSurf(surf,FC4h);
    %display(['FC4h = [' num2str(FC4h(1,1)) ' ' num2str(FC4h(1,2)) ' ' num2str(FC4h(1,3)) '], ' num2str(len) ' away from FT8']);
    
    [FC4,~] = curveWalk(curve_pts_FT7FCzFT8, 6*len_FT7FCzFT8/8);
    FC4 = bringPtsToSurf(surf,FC4);
    %display(['FC4 = [' num2str(FC4(1,1)) ' ' num2str(FC4(1,2)) ' ' num2str(FC4(1,3)) '], ' num2str(len) ' away from FT8']);
    
    [FC6h,~] = curveWalk(curve_pts_FT7FCzFT8, 6.5*len_FT7FCzFT8/8);
    FC6h = bringPtsToSurf(surf,FC6h);
    %display(['FC6h = [' num2str(FC6h(1,1)) ' ' num2str(FC6h(1,2)) ' ' num2str(FC6h(1,3)) '], ' num2str(len) ' away from FT8']);
    
    [FC6,~] = curveWalk(curve_pts_FT7FCzFT8, 7*len_FT7FCzFT8/8);
    FC6 = bringPtsToSurf(surf,FC6);
    %display(['FC6 = [' num2str(FC6(1,1)) ' ' num2str(FC6(1,2)) ' ' num2str(FC6(1,3)) '], ' num2str(len) ' away from FT8']);
    
    [FT8h,~] = curveWalk(curve_pts_FT7FCzFT8, 7.5*len_FT7FCzFT8/8);
    FT8h = bringPtsToSurf(surf,FT8h);
    %display(['FT8h = [' num2str(FT8h(1,1)) ' ' num2str(FT8h(1,2)) ' ' num2str(FT8h(1,3)) '], ' num2str(len) ' away from FT8']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find point along TP7-CPz-TP8 curve 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    plane = plane_equation(TP7, CPz, TP8);
    [curve_pts_TP7CPzTP8,len_TP7CPzTP8] = curveGen(TP7, TP8, CPz, plane, surf,dt);
    %display(['TP7-CPz-TP8 curve length: ' num2str(len_TP7CPzTP8)]);

    [TP7h,~] = curveWalk(curve_pts_TP7CPzTP8, 0.5*len_TP7CPzTP8/8);
    TP7h = bringPtsToSurf(surf,TP7h);
    %display(['TP7h = [' num2str(TP7h(1,1)) ' ' num2str(TP7h(1,2)) ' ' num2str(TP7h(1,3)) '], ' num2str(len) ' away from TP7']);
    
    [CP5,~] = curveWalk(curve_pts_TP7CPzTP8, 1*len_TP7CPzTP8/8);
    CP5 = bringPtsToSurf(surf,CP5);
    %display(['CP5 = [' num2str(CP5(1,1)) ' ' num2str(CP5(1,2)) ' ' num2str(CP5(1,3)) '], ' num2str(len) ' away from TP7']);
    
    [CP5h,~] = curveWalk(curve_pts_TP7CPzTP8, 1.5*len_TP7CPzTP8/8);
    CP5h = bringPtsToSurf(surf,CP5h);
    %display(['CP5h = [' num2str(CP5h(1,1)) ' ' num2str(CP5h(1,2)) ' ' num2str(CP5h(1,3)) '], ' num2str(len) ' away from TP7']);
    
    [CP3,~] = curveWalk(curve_pts_TP7CPzTP8, 2*len_TP7CPzTP8/8);
    CP3 = bringPtsToSurf(surf,CP3);
    %display(['CP3 = [' num2str(CP3(1,1)) ' ' num2str(CP3(1,2)) ' ' num2str(CP3(1,3)) '], ' num2str(len) ' away from TP7']);

    [CP3h,~] = curveWalk(curve_pts_TP7CPzTP8, 2.5*len_TP7CPzTP8/8);
    CP3h = bringPtsToSurf(surf,CP3h);
    %display(['CP3h = [' num2str(CP3h(1,1)) ' ' num2str(CP3h(1,2)) ' ' num2str(CP3h(1,3)) '], ' num2str(len) ' away from TP7']);
    
    [CP1,~] = curveWalk(curve_pts_TP7CPzTP8, 3*len_TP7CPzTP8/8);
    CP1 = bringPtsToSurf(surf,CP1);
    %display(['CP1 = [' num2str(CP1(1,1)) ' ' num2str(CP1(1,2)) ' ' num2str(CP1(1,3)) '], ' num2str(len) ' away from TP7']);
    
    [CP1h,~] = curveWalk(curve_pts_TP7CPzTP8, 3.5*len_TP7CPzTP8/8);
    CP1h = bringPtsToSurf(surf,CP1h);
    %display(['CP1h = [' num2str(CP1h(1,1)) ' ' num2str(CP1h(1,2)) ' ' num2str(CP1h(1,3)) '], ' num2str(len) ' away from TP7']);
    
    [CP2h,~] = curveWalk(curve_pts_TP7CPzTP8, 4.5*len_TP7CPzTP8/8);
    CP2h = bringPtsToSurf(surf,CP2h);
    %display(['CP2h = [' num2str(CP2h(1,1)) ' ' num2str(CP2h(1,2)) ' ' num2str(CP2h(1,3)) '], ' num2str(len) ' away from TP8']);
    
    [CP2,~] = curveWalk(curve_pts_TP7CPzTP8, 5*len_TP7CPzTP8/8);
    CP2 = bringPtsToSurf(surf,CP2);
    %display(['CP2 = [' num2str(CP2(1,1)) ' ' num2str(CP2(1,2)) ' ' num2str(CP2(1,3)) '], ' num2str(len) ' away from TP8']);
    
    [CP4h,~] = curveWalk(curve_pts_TP7CPzTP8, 5.5*len_TP7CPzTP8/8);
    CP4h = bringPtsToSurf(surf,CP4h);
    %display(['CP4h = [' num2str(CP4h(1,1)) ' ' num2str(CP4h(1,2)) ' ' num2str(CP4h(1,3)) '], ' num2str(len) ' away from TP8']);
    
    [CP4,~] = curveWalk(curve_pts_TP7CPzTP8, 6*len_TP7CPzTP8/8);
    CP4 = bringPtsToSurf(surf,CP4);
    %display(['CP4 = [' num2str(CP4(1,1)) ' ' num2str(CP4(1,2)) ' ' num2str(CP4(1,3)) '], ' num2str(len) ' away from TP8']);
    
    [CP6h,~] = curveWalk(curve_pts_TP7CPzTP8, 6.5*len_TP7CPzTP8/8);
    CP6h = bringPtsToSurf(surf,CP6h);
    %display(['CP6h = [' num2str(CP6h(1,1)) ' ' num2str(CP6h(1,2)) ' ' num2str(CP6h(1,3)) '], ' num2str(len) ' away from TP8']);
    
    [CP6,~] = curveWalk(curve_pts_TP7CPzTP8, 7*len_TP7CPzTP8/8);
    CP6 = bringPtsToSurf(surf,CP6);
    %display(['CP6 = [' num2str(CP6(1,1)) ' ' num2str(CP6(1,2)) ' ' num2str(CP6(1,3)) '], ' num2str(len) ' away from TP8']);
    
    [TP8h,~] = curveWalk(curve_pts_TP7CPzTP8, 7.5*len_TP7CPzTP8/8);
    TP8h = bringPtsToSurf(surf,TP8h);
    %display(['TP8h = [' num2str(TP8h(1,1)) ' ' num2str(TP8h(1,2)) ' ' num2str(TP8h(1,3)) '], ' num2str(len) ' away from TP8']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find point along FTT7-FCCz-FTT8 curve 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    plane = plane_equation(FTT7, FCCz, FTT8);
    [curve_pts_FTT7FCCzFTT8,len_FTT7FCCzFTT8] = curveGen(FTT7, FTT8, FCCz, plane, surf,dt);
    %display(['FTT7-FCCz-FTT8 curve length: ' num2str(len_FTT7FCCzFTT8)]);
    
    [FTT7h,~] = curveWalk(curve_pts_FTT7FCCzFTT8, 0.5*len_FTT7FCCzFTT8/8);
    FTT7h = bringPtsToSurf(surf,FTT7h);
    %display(['FTT7h = [' num2str(FTT7h(1,1)) ' ' num2str(FTT7h(1,2)) ' ' num2str(FTT7h(1,3)) '], ' num2str(len) ' away from FTT7']);

    [FCC5,~] = curveWalk(curve_pts_FTT7FCCzFTT8, 1*len_FTT7FCCzFTT8/8);
    FCC5 = bringPtsToSurf(surf,FCC5);
    %display(['FCC5 = [' num2str(FCC5(1,1)) ' ' num2str(FCC5(1,2)) ' ' num2str(FCC5(1,3)) '], ' num2str(len) ' away from FTT7']);
    
    [FCC5h,~] = curveWalk(curve_pts_FTT7FCCzFTT8, 1.5*len_FTT7FCCzFTT8/8);
    FCC5h = bringPtsToSurf(surf,FCC5h);
    %display(['FCC5h = [' num2str(FCC5h(1,1)) ' ' num2str(FCC5h(1,2)) ' ' num2str(FC5h(1,3)) '], ' num2str(len) ' away from FTT7']);
    
    [FCC3,~] = curveWalk(curve_pts_FTT7FCCzFTT8, 2*len_FTT7FCCzFTT8/8);
    FCC3 = bringPtsToSurf(surf,FCC3);
    %display(['FCC3 = [' num2str(FCC3(1,1)) ' ' num2str(FCC3(1,2)) ' ' num2str(FCC3(1,3)) '], ' num2str(len) ' away from FTT7']);

    [FCC3h,~] = curveWalk(curve_pts_FTT7FCCzFTT8, 2.5*len_FTT7FCCzFTT8/8);
    FCC3h = bringPtsToSurf(surf,FCC3h);
    %display(['FC3h = [' num2str(FCC3h(1,1)) ' ' num2str(FCC3h(1,2)) ' ' num2str(FCC3h(1,3)) '], ' num2str(len) ' away from FTT7']);
    
    [FCC1,~] = curveWalk(curve_pts_FTT7FCCzFTT8, 3*len_FTT7FCCzFTT8/8);
    FCC1 = bringPtsToSurf(surf,FCC1);
    %display(['FCC1 = [' num2str(FCC1(1,1)) ' ' num2str(FCC1(1,2)) ' ' num2str(FCC1(1,3)) '], ' num2str(len) ' away from FTT7']);
    
    [FCC1h,~] = curveWalk(curve_pts_FTT7FCCzFTT8, 3.5*len_FTT7FCCzFTT8/8);
    FCC1h = bringPtsToSurf(surf,FCC1h);
    %display(['FCC1h = [' num2str(FCC1h(1,1)) ' ' num2str(FCC1h(1,2)) ' ' num2str(FCC1h(1,3)) '], ' num2str(len) ' away from FTT7']);
    
    [FCC2h,~] = curveWalk(curve_pts_FTT7FCCzFTT8, 4.5*len_FTT7FCCzFTT8/8);
    FCC2h = bringPtsToSurf(surf,FCC2h);
    %display(['FCC2h = [' num2str(FCC2h(1,1)) ' ' num2str(FCC2h(1,2)) ' ' num2str(FCC2h(1,3)) '], ' num2str(len) ' away from FTT8']);
    
    [FCC2,~] = curveWalk(curve_pts_FTT7FCCzFTT8, 5*len_FTT7FCCzFTT8/8);
    FCC2 = bringPtsToSurf(surf,FCC2);
    %display(['FCC2 = [' num2str(FCC2(1,1)) ' ' num2str(FCC2(1,2)) ' ' num2str(FCC2(1,3)) '], ' num2str(len) ' away from FTT8']);
    
    [FCC4h,~] = curveWalk(curve_pts_FTT7FCCzFTT8, 5.5*len_FTT7FCCzFTT8/8);
    FCC4h = bringPtsToSurf(surf,FCC4h);
    %display(['FCC4h = [' num2str(FCC4h(1,1)) ' ' num2str(FCC4h(1,2)) ' ' num2str(FCC4h(1,3)) '], ' num2str(len) ' away from FTT8']);
    
    [FCC4,~] = curveWalk(curve_pts_FTT7FCCzFTT8, 6*len_FTT7FCCzFTT8/8);
    FCC4 = bringPtsToSurf(surf,FCC4);
    %display(['FCC4 = [' num2str(FCC4(1,1)) ' ' num2str(FCC4(1,2)) ' ' num2str(FCC4(1,3)) '], ' num2str(len) ' away from FTT8']);
    
    [FCC6h,~] = curveWalk(curve_pts_FTT7FCCzFTT8, 6.5*len_FTT7FCCzFTT8/8);
    FCC6h = bringPtsToSurf(surf,FCC6h);
    %display(['FCC6h = [' num2str(FCC6h(1,1)) ' ' num2str(FCC6h(1,2)) ' ' num2str(FCC6h(1,3)) '], ' num2str(len) ' away from FTT8']);
    
    [FCC6,~] = curveWalk(curve_pts_FTT7FCCzFTT8, 7*len_FTT7FCCzFTT8/8);
    FCC6 = bringPtsToSurf(surf,FCC6);
    %display(['FCC6 = [' num2str(FCC6(1,1)) ' ' num2str(FCC6(1,2)) ' ' num2str(FCC6(1,3)) '], ' num2str(len) ' away from FTT8']);
    
    [FTT8h,~] = curveWalk(curve_pts_FTT7FCCzFTT8, 7.5*len_FTT7FCCzFTT8/8);
    FTT8h = bringPtsToSurf(surf,FTT8h);
    %display(['FTT8h = [' num2str(FTT8h(1,1)) ' ' num2str(FTT8h(1,2)) ' ' num2str(FTT8h(1,3)) '], ' num2str(len) ' away from FTT8']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find point along TTP7-CCPz-TTP8 curve 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    plane = plane_equation(TTP7, CCPz, TTP8);
    [curve_pts_TTP7CCPzTTP8,len_TTP7CCPzTTP8] = curveGen(TTP7, TTP8, CCPz, plane, surf,dt);
    %display(['TTP7-CCPz-TTP8 curve length: ' num2str(len_TTP7CCPzTTP8)]);

    [TTP7h,~] = curveWalk(curve_pts_TTP7CCPzTTP8, 0.5*len_TTP7CCPzTTP8/8);
    TTP7h = bringPtsToSurf(surf,TTP7h);
    %display(['TTP7h = [' num2str(TTP7h(1,1)) ' ' num2str(TTP7h(1,2)) ' ' num2str(TTP7h(1,3)) '], ' num2str(len) ' away from TTP7']);
    
    [CCP5,~] = curveWalk(curve_pts_TTP7CCPzTTP8, 1*len_TTP7CCPzTTP8/8);
    CCP5 = bringPtsToSurf(surf,CCP5);
    %display(['CCP5 = [' num2str(CCP5(1,1)) ' ' num2str(CCP5(1,2)) ' ' num2str(CCP5(1,3)) '], ' num2str(len) ' away from TTP7']);
    
    [CCP5h,~] = curveWalk(curve_pts_TTP7CCPzTTP8, 1.5*len_TTP7CCPzTTP8/8);
    CCP5h = bringPtsToSurf(surf,CCP5h);
    %display(['CCP5h = [' num2str(CCP5h(1,1)) ' ' num2str(CCP5h(1,2)) ' ' num2str(CCP5h(1,3)) '], ' num2str(len) ' away from TTP7']);
    
    [CCP3,~] = curveWalk(curve_pts_TTP7CCPzTTP8, 2*len_TTP7CCPzTTP8/8);
    CCP3 = bringPtsToSurf(surf,CCP3);
    %display(['CCP3 = [' num2str(CCP3(1,1)) ' ' num2str(CCP3(1,2)) ' ' num2str(CCP3(1,3)) '], ' num2str(len) ' away from TTP7']);

    [CCP3h,~] = curveWalk(curve_pts_TTP7CCPzTTP8, 2.5*len_TTP7CCPzTTP8/8);
    CCP3h = bringPtsToSurf(surf,CCP3h);
    %display(['CCP3h = [' num2str(CCP3h(1,1)) ' ' num2str(CCP3h(1,2)) ' ' num2str(CCP3h(1,3)) '], ' num2str(len) ' away from TTP7']);
    
    [CCP1,~] = curveWalk(curve_pts_TTP7CCPzTTP8, 3*len_TTP7CCPzTTP8/8);
    CCP1 = bringPtsToSurf(surf,CCP1);
    %display(['CP1 = [' num2str(CCP1(1,1)) ' ' num2str(CCP1(1,2)) ' ' num2str(CCP1(1,3)) '], ' num2str(len) ' away from TTP7']);
    
    [CCP1h,~] = curveWalk(curve_pts_TTP7CCPzTTP8, 3.5*len_TTP7CCPzTTP8/8);
    CCP1h = bringPtsToSurf(surf,CCP1h);
    %display(['CCP1h = [' num2str(CCP1h(1,1)) ' ' num2str(CCP1h(1,2)) ' ' num2str(CCP1h(1,3)) '], ' num2str(len) ' away from TTP7']);
    
    [CCP2h,~] = curveWalk(curve_pts_TTP7CCPzTTP8, 4.5*len_TTP7CCPzTTP8/8);
    CCP2h = bringPtsToSurf(surf,CCP2h);
    %display(['CCP2h = [' num2str(CCP2h(1,1)) ' ' num2str(CCP2h(1,2)) ' ' num2str(CCP2h(1,3)) '], ' num2str(len) ' away from TTP8']);
    
    [CCP2,~] = curveWalk(curve_pts_TTP7CCPzTTP8, 5*len_TTP7CCPzTTP8/8);
    CCP2 = bringPtsToSurf(surf,CCP2);
    %display(['CCP2 = [' num2str(CCP2(1,1)) ' ' num2str(CCP2(1,2)) ' ' num2str(CCP2(1,3)) '], ' num2str(len) ' away from TTP8']);
    
    [CCP4h,~] = curveWalk(curve_pts_TTP7CCPzTTP8, 5.5*len_TTP7CCPzTTP8/8);
    CCP4h = bringPtsToSurf(surf,CCP4h);
    %display(['CCP4h = [' num2str(CCP4h(1,1)) ' ' num2str(CCP4h(1,2)) ' ' num2str(CCP4h(1,3)) '], ' num2str(len) ' away from TTP8']);
    
    [CCP4,~] = curveWalk(curve_pts_TTP7CCPzTTP8, 6*len_TTP7CCPzTTP8/8);
    CCP4 = bringPtsToSurf(surf,CCP4);
    %display(['CCP4 = [' num2str(CCP4(1,1)) ' ' num2str(CCP4(1,2)) ' ' num2str(CCP4(1,3)) '], ' num2str(len) ' away from TTP8']);
    
    [CCP6h,~] = curveWalk(curve_pts_TTP7CCPzTTP8, 6.5*len_TTP7CCPzTTP8/8);
    CCP6h = bringPtsToSurf(surf,CCP6h);
    %display(['CCP6h = [' num2str(CCP6h(1,1)) ' ' num2str(CCP6h(1,2)) ' ' num2str(CCP6h(1,3)) '], ' num2str(len) ' away from TTP8']);
    
    [CCP6,~] = curveWalk(curve_pts_TTP7CCPzTTP8, 7*len_TTP7CCPzTTP8/8);
    CCP6 = bringPtsToSurf(surf,CCP6);
    %display(['CCP6 = [' num2str(CCP6(1,1)) ' ' num2str(CCP6(1,2)) ' ' num2str(CCP6(1,3)) '], ' num2str(len) ' away from TTP8']);
    
    [TTP8h,~] = curveWalk(curve_pts_TTP7CCPzTTP8, 7.5*len_TTP7CCPzTTP8/8);
    TTP8h = bringPtsToSurf(surf,TTP8h);
    %display(['TTP8h = [' num2str(TTP8h(1,1)) ' ' num2str(TTP8h(1,2)) ' ' num2str(TTP8h(1,3)) '], ' num2str(len) ' away from TTP8']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find point along FFT7-FFCz-FFT8 curve 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    plane = plane_equation(FFT7, FFCz, FFT8);
    [curve_pts_FFT7FFCzFFT8,len_FFT7FFCzFFT8] = curveGen(FFT7, FFT8, FFCz, plane, surf,dt);
    %display(['FFT7-FFCz-FFT8 curve length: ' num2str(len_FFT7FFCzFFT8)]);
    
    [FFT7h,~] = curveWalk(curve_pts_FFT7FFCzFFT8, 0.5*len_FFT7FFCzFFT8/8);
    FFT7h = bringPtsToSurf(surf,FFT7h);
    %display(['FFT7h = [' num2str(FFT7h(1,1)) ' ' num2str(FFT7h(1,2)) ' ' num2str(FFT7h(1,3)) '], ' num2str(len) ' away from FFT7']);

    [FFC5,~] = curveWalk(curve_pts_FFT7FFCzFFT8, 1*len_FFT7FFCzFFT8/8);
    FFC5 = bringPtsToSurf(surf,FFC5);
    %display(['FFC5 = [' num2str(FFC5(1,1)) ' ' num2str(FFC5(1,2)) ' ' num2str(FFC5(1,3)) '], ' num2str(len) ' away from FFT7']);
    
    [FFC5h,~] = curveWalk(curve_pts_FFT7FFCzFFT8, 1.5*len_FFT7FFCzFFT8/8);
    FFC5h = bringPtsToSurf(surf,FFC5h);
    %display(['FFC5h = [' num2str(FFC5h(1,1)) ' ' num2str(FFC5h(1,2)) ' ' num2str(FFC5h(1,3)) '], ' num2str(len) ' away from FFT7']);
    
    [FFC3,~] = curveWalk(curve_pts_FFT7FFCzFFT8, 2*len_FFT7FFCzFFT8/8);
    FFC3 = bringPtsToSurf(surf,FFC3);
    %display(['FFC3 = [' num2str(FFC3(1,1)) ' ' num2str(FFC3(1,2)) ' ' num2str(FFC3(1,3)) '], ' num2str(len) ' away from FFT7']);

    [FFC3h,~] = curveWalk(curve_pts_FFT7FFCzFFT8, 2.5*len_FFT7FFCzFFT8/8);
    FFC3h = bringPtsToSurf(surf,FFC3h);
    %display(['FFC3h = [' num2str(FFC3h(1,1)) ' ' num2str(FFC3h(1,2)) ' ' num2str(FFC3h(1,3)) '], ' num2str(len) ' away from FFT7']);
    
    [FFC1,~] = curveWalk(curve_pts_FFT7FFCzFFT8, 3*len_FFT7FFCzFFT8/8);
    FFC1 = bringPtsToSurf(surf,FFC1);
    %display(['FFC1 = [' num2str(FFC1(1,1)) ' ' num2str(FFC1(1,2)) ' ' num2str(FFC1(1,3)) '], ' num2str(len) ' away from FFT7']);
    
    [FFC1h,~] = curveWalk(curve_pts_FFT7FFCzFFT8, 3.5*len_FFT7FFCzFFT8/8);
    FFC1h = bringPtsToSurf(surf,FFC1h);
    %display(['FFC1h = [' num2str(FFC1h(1,1)) ' ' num2str(FFC1h(1,2)) ' ' num2str(FFC1h(1,3)) '], ' num2str(len) ' away from FFT7']);
    
    [FFC2h,~] = curveWalk(curve_pts_FFT7FFCzFFT8, 4.5*len_FFT7FFCzFFT8/8);
    FFC2h = bringPtsToSurf(surf,FFC2h);
    %display(['FFC2h = [' num2str(FFC2h(1,1)) ' ' num2str(FFC2h(1,2)) ' ' num2str(FFC2h(1,3)) '], ' num2str(len) ' away from FFT8']);
    
    [FFC2,~] = curveWalk(curve_pts_FFT7FFCzFFT8, 5*len_FFT7FFCzFFT8/8);
    FFC2 = bringPtsToSurf(surf,FFC2);
    %display(['FFC2 = [' num2str(FFC2(1,1)) ' ' num2str(FFC2(1,2)) ' ' num2str(FFC2(1,3)) '], ' num2str(len) ' away from FFT8']);
    
    [FFC4h,~] = curveWalk(curve_pts_FFT7FFCzFFT8, 5.5*len_FFT7FFCzFFT8/8);
    FFC4h = bringPtsToSurf(surf,FFC4h);
    %display(['FFC4h = [' num2str(FFC4h(1,1)) ' ' num2str(FFC4h(1,2)) ' ' num2str(FFC4h(1,3)) '], ' num2str(len) ' away from FFT8']);
    
    [FFC4,~] = curveWalk(curve_pts_FFT7FFCzFFT8, 6*len_FFT7FFCzFFT8/8);
    FFC4 = bringPtsToSurf(surf,FFC4);
    %display(['FFC4 = [' num2str(FFC4(1,1)) ' ' num2str(FFC4(1,2)) ' ' num2str(FFC4(1,3)) '], ' num2str(len) ' away from FFT8']);
    
    [FFC6h,~] = curveWalk(curve_pts_FFT7FFCzFFT8, 6.5*len_FFT7FFCzFFT8/8);
    FFC6h = bringPtsToSurf(surf,FFC6h);
    %display(['FFC6h = [' num2str(FFC6h(1,1)) ' ' num2str(FFC6h(1,2)) ' ' num2str(FFC6h(1,3)) '], ' num2str(len) ' away from FFT8']);
    
    [FFC6,~] = curveWalk(curve_pts_FFT7FFCzFFT8, 7*len_FFT7FFCzFFT8/8);
    FFC6 = bringPtsToSurf(surf,FFC6);
    %display(['FFC6 = [' num2str(FFC6(1,1)) ' ' num2str(FFC6(1,2)) ' ' num2str(FFC6(1,3)) '], ' num2str(len) ' away from FFT8']);
    
    [FFT8h,~] = curveWalk(curve_pts_FFT7FFCzFFT8, 7.5*len_FFT7FFCzFFT8/8);
    FFT8h = bringPtsToSurf(surf,FFT8h);
    %display(['FFT8h = [' num2str(FFT8h(1,1)) ' ' num2str(FFT8h(1,2)) ' ' num2str(FFT8h(1,3)) '], ' num2str(len) ' away from FFT8']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find point along TPP7-CPPz-TPP8 curve 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    plane = plane_equation(TPP7, CPPz, TPP8);
    [curve_pts_TPP7CPPzTPP8,len_TPP7CPPzTPP8] = curveGen(TPP7, TPP8, CPPz, plane, surf,dt);
    %display(['TPP7-CPPz-TPP8 curve length: ' num2str(len_TPP7CPPzTPP8)]);

    [TPP7h,~] = curveWalk(curve_pts_TPP7CPPzTPP8, 0.5*len_TPP7CPPzTPP8/8);
    TPP7h = bringPtsToSurf(surf,TPP7h);
    %display(['TPP7h = [' num2str(TPP7h(1,1)) ' ' num2str(TPP7h(1,2)) ' ' num2str(TPP7h(1,3)) '], ' num2str(len) ' away from TPP7']);
    
    [CPP5,~] = curveWalk(curve_pts_TPP7CPPzTPP8, 1*len_TPP7CPPzTPP8/8);
    CPP5 = bringPtsToSurf(surf,CPP5);
    %display(['CPP5 = [' num2str(CPP5(1,1)) ' ' num2str(CPP5(1,2)) ' ' num2str(CPP5(1,3)) '], ' num2str(len) ' away from TPP7']);
    
    [CPP5h,~] = curveWalk(curve_pts_TPP7CPPzTPP8, 1.5*len_TPP7CPPzTPP8/8);
    CPP5h = bringPtsToSurf(surf,CPP5h);
    %display(['CPP5h = [' num2str(CPP5h(1,1)) ' ' num2str(CPP5h(1,2)) ' ' num2str(CPP5h(1,3)) '], ' num2str(len) ' away from TPP7']);
    
    [CPP3,~] = curveWalk(curve_pts_TPP7CPPzTPP8, 2*len_TPP7CPPzTPP8/8);
    CPP3 = bringPtsToSurf(surf,CPP3);
    %display(['CPP3 = [' num2str(CPP3(1,1)) ' ' num2str(CPP3(1,2)) ' ' num2str(CPP3(1,3)) '], ' num2str(len) ' away from TPP7']);

    [CPP3h,~] = curveWalk(curve_pts_TPP7CPPzTPP8, 2.5*len_TPP7CPPzTPP8/8);
    CPP3h = bringPtsToSurf(surf,CPP3h);
    %display(['CPP3h = [' num2str(CPP3h(1,1)) ' ' num2str(CPP3h(1,2)) ' ' num2str(CPP3h(1,3)) '], ' num2str(len) ' away from TPP7']);
    
    [CPP1,~] = curveWalk(curve_pts_TPP7CPPzTPP8, 3*len_TPP7CPPzTPP8/8);
    CPP1 = bringPtsToSurf(surf,CPP1);
    %display(['CPP1 = [' num2str(CPP1(1,1)) ' ' num2str(CPP1(1,2)) ' ' num2str(CPP1(1,3)) '], ' num2str(len) ' away from TPP7']);
    
    [CPP1h,~] = curveWalk(curve_pts_TPP7CPPzTPP8, 3.5*len_TPP7CPPzTPP8/8);
    CPP1h = bringPtsToSurf(surf,CPP1h);
    %display(['CPP1h = [' num2str(CPP1h(1,1)) ' ' num2str(CPP1h(1,2)) ' ' num2str(CPP1h(1,3)) '], ' num2str(len) ' away from TPP7']);
    
    [CPP2h,~] = curveWalk(curve_pts_TPP7CPPzTPP8, 4.5*len_TPP7CPPzTPP8/8);
    CPP2h = bringPtsToSurf(surf,CPP2h);
    %display(['CPP2h = [' num2str(CPP2h(1,1)) ' ' num2str(CPP2h(1,2)) ' ' num2str(CPP2h(1,3)) '], ' num2str(len) ' away from TPP8']);
    
    [CPP2,~] = curveWalk(curve_pts_TPP7CPPzTPP8, 5*len_TPP7CPPzTPP8/8);
    CPP2 = bringPtsToSurf(surf,CPP2);
    %display(['CPP2 = [' num2str(CPP2(1,1)) ' ' num2str(CPP2(1,2)) ' ' num2str(CPP2(1,3)) '], ' num2str(len) ' away from TPP8']);
    
    [CPP4h,~] = curveWalk(curve_pts_TPP7CPPzTPP8, 5.5*len_TPP7CPPzTPP8/8);
    CPP4h = bringPtsToSurf(surf,CPP4h);
    %display(['CPP4h = [' num2str(CPP4h(1,1)) ' ' num2str(CPP4h(1,2)) ' ' num2str(CPP4h(1,3)) '], ' num2str(len) ' away from TPP8']);
    
    [CPP4,~] = curveWalk(curve_pts_TPP7CPPzTPP8, 6*len_TPP7CPPzTPP8/8);
    CPP4 = bringPtsToSurf(surf,CPP4);
    %display(['CPP4 = [' num2str(CPP4(1,1)) ' ' num2str(CPP4(1,2)) ' ' num2str(CPP4(1,3)) '], ' num2str(len) ' away from TPP8']);
    
    [CPP6h,~] = curveWalk(curve_pts_TPP7CPPzTPP8, 6.5*len_TPP7CPPzTPP8/8);
    CPP6h = bringPtsToSurf(surf,CPP6h);
    %display(['CPP6h = [' num2str(CPP6h(1,1)) ' ' num2str(CPP6h(1,2)) ' ' num2str(CPP6h(1,3)) '], ' num2str(len) ' away from TPP8']);
    
    [CPP6,~] = curveWalk(curve_pts_TPP7CPPzTPP8, 7*len_TPP7CPPzTPP8/8);
    CPP6 = bringPtsToSurf(surf,CPP6);
    %display(['CPP6 = [' num2str(CPP6(1,1)) ' ' num2str(CPP6(1,2)) ' ' num2str(CPP6(1,3)) '], ' num2str(len) ' away from TPP8']);
    
    [TPP8h,~] = curveWalk(curve_pts_TPP7CPPzTPP8, 7.5*len_TPP7CPPzTPP8/8);
    TPP8h = bringPtsToSurf(surf,TPP8h);
    %display(['TPP8h = [' num2str(TPP8h(1,1)) ' ' num2str(TPP8h(1,2)) ' ' num2str(TPP8h(1,3)) '], ' num2str(len) ' away from TPP8']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find point along AFF7-AFFz-AFF8 curve 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    plane = plane_equation(AFF7, AFFz, AFF8);
    [curve_pts_AFF7AFFzAFF8,len_AFF7AFFzAFF8] = curveGen(AFF7, AFF8, AFFz, plane, surf,dt);
    %display(['AFF7-AFFz-AFF8 curve length: ' num2str(len_AFF7AFFzAFF8)]);
    
    [AFF7h,~] = curveWalk(curve_pts_AFF7AFFzAFF8, 0.5*len_AFF7AFFzAFF8/8);
    AFF7h = bringPtsToSurf(surf,AFF7h);
    %display(['AFF7h = [' num2str(AFF7h(1,1)) ' ' num2str(AFF7h(1,2)) ' ' num2str(AFF7h(1,3)) '], ' num2str(len) ' away from AFF7']);

    [AFF5,~] = curveWalk(curve_pts_AFF7AFFzAFF8, 1*len_AFF7AFFzAFF8/8);
    AFF5 = bringPtsToSurf(surf,AFF5);
    %display(['AFF5 = [' num2str(AFF5(1,1)) ' ' num2str(AFF5(1,2)) ' ' num2str(AFF5(1,3)) '], ' num2str(len) ' away from AFF7']);
    
    [AFF5h,~] = curveWalk(curve_pts_AFF7AFFzAFF8, 1.5*len_AFF7AFFzAFF8/8);
    AFF5h = bringPtsToSurf(surf,AFF5h);
    %display(['FFC5h = [' num2str(AFF5h(1,1)) ' ' num2str(AFF5h(1,2)) ' ' num2str(AFF5h(1,3)) '], ' num2str(len) ' away from AFF7']);
    
    [AFF3,~] = curveWalk(curve_pts_AFF7AFFzAFF8, 2*len_AFF7AFFzAFF8/8);
    AFF3 = bringPtsToSurf(surf,AFF3);
    %display(['AFF3 = [' num2str(AFF3(1,1)) ' ' num2str(AFF3(1,2)) ' ' num2str(AFF3(1,3)) '], ' num2str(len) ' away from AFF7']);

    [AFF3h,~] = curveWalk(curve_pts_AFF7AFFzAFF8, 2.5*len_AFF7AFFzAFF8/8);
    AFF3h = bringPtsToSurf(surf,AFF3h);
    %display(['AFF3h = [' num2str(AFF3h(1,1)) ' ' num2str(AFF3h(1,2)) ' ' num2str(AFF3h(1,3)) '], ' num2str(len) ' away from AFF7']);
    
    [AFF1,~] = curveWalk(curve_pts_AFF7AFFzAFF8, 3*len_AFF7AFFzAFF8/8);
    AFF1 = bringPtsToSurf(surf,AFF1);
    %display(['AFF1 = [' num2str(AFF1(1,1)) ' ' num2str(AFF1(1,2)) ' ' num2str(AFF1(1,3)) '], ' num2str(len) ' away from AFF7']);
    
    [AFF1h,~] = curveWalk(curve_pts_AFF7AFFzAFF8, 3.5*len_AFF7AFFzAFF8/8);
    AFF1h = bringPtsToSurf(surf,AFF1h);
    %display(['AFF1h = [' num2str(AFF1h(1,1)) ' ' num2str(AFF1h(1,2)) ' ' num2str(AFF1h(1,3)) '], ' num2str(len) ' away from AFF7']);
    
    [AFF2h,~] = curveWalk(curve_pts_AFF7AFFzAFF8, 4.5*len_AFF7AFFzAFF8/8);
    AFF2h = bringPtsToSurf(surf,AFF2h);
    %display(['AFF2h = [' num2str(AFF2h(1,1)) ' ' num2str(AFF2h(1,2)) ' ' num2str(AFF2h(1,3)) '], ' num2str(len) ' away from AFF8']);
    
    [AFF2,~] = curveWalk(curve_pts_AFF7AFFzAFF8, 5*len_AFF7AFFzAFF8/8);
    AFF2 = bringPtsToSurf(surf,AFF2);
    %display(['AFF2 = [' num2str(AFF2(1,1)) ' ' num2str(AFF2(1,2)) ' ' num2str(AFF2(1,3)) '], ' num2str(len) ' away from AFF8']);
    
    [AFF4h,~] = curveWalk(curve_pts_AFF7AFFzAFF8, 5.5*len_AFF7AFFzAFF8/8);
    AFF4h = bringPtsToSurf(surf,AFF4h);
    %display(['AFF4h = [' num2str(AFF4h(1,1)) ' ' num2str(AFF4h(1,2)) ' ' num2str(AFF4h(1,3)) '], ' num2str(len) ' away from AFF8']);
    
    [AFF4,~] = curveWalk(curve_pts_AFF7AFFzAFF8, 6*len_AFF7AFFzAFF8/8);
    AFF4 = bringPtsToSurf(surf,AFF4);
    %display(['AFF4 = [' num2str(AFF4(1,1)) ' ' num2str(AFF4(1,2)) ' ' num2str(AFF4(1,3)) '], ' num2str(len) ' away from AFF8']);
    
    [AFF6h,~] = curveWalk(curve_pts_AFF7AFFzAFF8, 6.5*len_AFF7AFFzAFF8/8);
    AFF6h = bringPtsToSurf(surf,AFF6h);
    %display(['AFF6h = [' num2str(AFF6h(1,1)) ' ' num2str(AFF6h(1,2)) ' ' num2str(AFF6h(1,3)) '], ' num2str(len) ' away from AFF8']);
    
    [AFF6,~] = curveWalk(curve_pts_AFF7AFFzAFF8, 7*len_AFF7AFFzAFF8/8);
    AFF6 = bringPtsToSurf(surf,AFF6);
    %display(['AFF6 = [' num2str(AFF6(1,1)) ' ' num2str(AFF6(1,2)) ' ' num2str(AFF6(1,3)) '], ' num2str(len) ' away from AFF8']);
    
    [AFF8h,~] = curveWalk(curve_pts_AFF7AFFzAFF8, 7.5*len_AFF7AFFzAFF8/8);
    AFF8h = bringPtsToSurf(surf,AFF8h);
    %display(['AFF8h = [' num2str(AFF8h(1,1)) ' ' num2str(AFF8h(1,2)) ' ' num2str(AFF8h(1,3)) '], ' num2str(len) ' away from AFF8']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find point along PPO7-PPOz-PPO8 curve 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    plane = plane_equation(PPO7, PPOz, PPO8);
    [curve_pts_PPO7PPOzPPO8,len_PPO7PPOzPPO8] = curveGen(PPO7, PPO8, PPOz, plane, surf,dt);
    %display(['PPO7-PPOz-PPO8 curve length: ' num2str(len_PPO7PPOzPPO8)]);

    [PPO7h,~] = curveWalk(curve_pts_PPO7PPOzPPO8, 0.5*len_PPO7PPOzPPO8/8);
    PPO7h = bringPtsToSurf(surf,PPO7h);
    %display(['PPO7h = [' num2str(PPO7h(1,1)) ' ' num2str(PPO7h(1,2)) ' ' num2str(PPO7h(1,3)) '], ' num2str(len) ' away from PPO7']);
    
    [PPO5,~] = curveWalk(curve_pts_PPO7PPOzPPO8, 1*len_PPO7PPOzPPO8/8);
    PPO5 = bringPtsToSurf(surf,PPO5);
    %display(['PPO5 = [' num2str(PPO5(1,1)) ' ' num2str(PPO5(1,2)) ' ' num2str(PPO5(1,3)) '], ' num2str(len) ' away from PPO7']);
    
    [PPO5h,~] = curveWalk(curve_pts_PPO7PPOzPPO8, 1.5*len_PPO7PPOzPPO8/8);
    PPO5h = bringPtsToSurf(surf,PPO5h);
    %display(['PPO5h = [' num2str(PPO5h(1,1)) ' ' num2str(PPO5h(1,2)) ' ' num2str(PPO5h(1,3)) '], ' num2str(len) ' away from PPO7']);
    
    [PPO3,~] = curveWalk(curve_pts_PPO7PPOzPPO8, 2*len_PPO7PPOzPPO8/8);
    PPO3 = bringPtsToSurf(surf,PPO3);
    %display(['PPO3 = [' num2str(PPO3(1,1)) ' ' num2str(PPO3(1,2)) ' ' num2str(PPO3(1,3)) '], ' num2str(len) ' away from PPO7']);

    [PPO3h,~] = curveWalk(curve_pts_PPO7PPOzPPO8, 2.5*len_PPO7PPOzPPO8/8);
    PPO3h = bringPtsToSurf(surf,PPO3h);
    %display(['PPO3h = [' num2str(PPO3h(1,1)) ' ' num2str(PPO3h(1,2)) ' ' num2str(PPO3h(1,3)) '], ' num2str(len) ' away from PPO7']);
    
    [PPO1,~] = curveWalk(curve_pts_PPO7PPOzPPO8, 3*len_PPO7PPOzPPO8/8);
    PPO1 = bringPtsToSurf(surf,PPO1);
    %display(['PPO1 = [' num2str(PPO1(1,1)) ' ' num2str(PPO1(1,2)) ' ' num2str(PPO1(1,3)) '], ' num2str(len) ' away from PPO7']);
    
    [PPO1h,~] = curveWalk(curve_pts_PPO7PPOzPPO8, 3.5*len_PPO7PPOzPPO8/8);
    PPO1h = bringPtsToSurf(surf,PPO1h);
    %display(['PPO1h = [' num2str(PPO1h(1,1)) ' ' num2str(PPO1h(1,2)) ' ' num2str(PPO1h(1,3)) '], ' num2str(len) ' away from PPO7']);
    
    [PPO2h,~] = curveWalk(curve_pts_PPO7PPOzPPO8, 4.5*len_PPO7PPOzPPO8/8);
    PPO2h = bringPtsToSurf(surf,PPO2h);
    %display(['PPO2h = [' num2str(PPO2h(1,1)) ' ' num2str(PPO2h(1,2)) ' ' num2str(PPO2h(1,3)) '], ' num2str(len) ' away from PPO8']);
    
    [PPO2,~] = curveWalk(curve_pts_PPO7PPOzPPO8, 5*len_PPO7PPOzPPO8/8);
    PPO2 = bringPtsToSurf(surf,PPO2);
    %display(['PPO2 = [' num2str(PPO2(1,1)) ' ' num2str(PPO2(1,2)) ' ' num2str(PPO2(1,3)) '], ' num2str(len) ' away from PPO8']);
    
    [PPO4h,~] = curveWalk(curve_pts_PPO7PPOzPPO8, 5.5*len_PPO7PPOzPPO8/8);
    PPO4h = bringPtsToSurf(surf,PPO4h);
    %display(['PPO4h = [' num2str(PPO4h(1,1)) ' ' num2str(PPO4h(1,2)) ' ' num2str(PPO4h(1,3)) '], ' num2str(len) ' away from PPO8']);
    
    [PPO4,~] = curveWalk(curve_pts_PPO7PPOzPPO8, 6*len_PPO7PPOzPPO8/8);
    PPO4 = bringPtsToSurf(surf,PPO4);
    %display(['PPO4 = [' num2str(PPO4(1,1)) ' ' num2str(PPO4(1,2)) ' ' num2str(PPO4(1,3)) '], ' num2str(len) ' away from PPO8']);
    
    [PPO6h,~] = curveWalk(curve_pts_PPO7PPOzPPO8, 6.5*len_PPO7PPOzPPO8/8);
    PPO6h = bringPtsToSurf(surf,PPO6h);
    %display(['PPO6h = [' num2str(PPO6h(1,1)) ' ' num2str(PPO6h(1,2)) ' ' num2str(PPO6h(1,3)) '], ' num2str(len) ' away from PPO8']);
    
    [PPO6,~] = curveWalk(curve_pts_PPO7PPOzPPO8, 7*len_PPO7PPOzPPO8/8);
    PPO6 = bringPtsToSurf(surf,PPO6);
    %display(['PPO6 = [' num2str(PPO6(1,1)) ' ' num2str(PPO6(1,2)) ' ' num2str(PPO6(1,3)) '], ' num2str(len) ' away from PPO8']);
    
    [PPO8h,~] = curveWalk(curve_pts_PPO7PPOzPPO8, 7.5*len_PPO7PPOzPPO8/8);
    PPO8h = bringPtsToSurf(surf,PPO8h);
    %display(['PPO8h = [' num2str(PPO8h(1,1)) ' ' num2str(PPO8h(1,2)) ' ' num2str(PPO8h(1,3)) '], ' num2str(len) ' away from PPO8']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find point along AFp7-AFpz-AFp8 curve 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    plane = plane_equation(AFp7, AFpz, AFp8);
    [curve_pts_AFp7AFpzAFp8,len_AFp7AFpzAFp8] = curveGen(AFp7, AFp8, AFpz, plane, surf,dt);
    %display(['AFp7-AFpz-AFp8 curve length: ' num2str(len_AFp7AFpzAFp8)]);

    [AFp5,~] = curveWalk(curve_pts_AFp7AFpzAFp8, 1*len_AFp7AFpzAFp8/8);
    AFp5 = bringPtsToSurf(surf,AFp5);
    %display(['AFp5 = [' num2str(AFp5(1,1)) ' ' num2str(AFp5(1,2)) ' ' num2str(AFp5(1,3)) '], ' num2str(len) ' away from AFp7']);
    
    [AFp3,~] = curveWalk(curve_pts_AFp7AFpzAFp8, 2*len_AFp7AFpzAFp8/8);
    AFp3 = bringPtsToSurf(surf,AFp3);
    %display(['AFp3 = [' num2str(AFp3(1,1)) ' ' num2str(AFp3(1,2)) ' ' num2str(AFp3(1,3)) '], ' num2str(len) ' away from AFp7']);

    [AFp1,~] = curveWalk(curve_pts_AFp7AFpzAFp8, 3*len_AFp7AFpzAFp8/8);
    AFp1 = bringPtsToSurf(surf,AFp1);
    %display(['AFp1 = [' num2str(AFp1(1,1)) ' ' num2str(AFp1(1,2)) ' ' num2str(AFp1(1,3)) '], ' num2str(len) ' away from AFp7']);
    
    [AFp2,~] = curveWalk(curve_pts_AFp7AFpzAFp8, 5*len_AFp7AFpzAFp8/8);
    AFp2 = bringPtsToSurf(surf,AFp2);
    %display(['AFp2 = [' num2str(AFp2(1,1)) ' ' num2str(AFp2(1,2)) ' ' num2str(AFp2(1,3)) '], ' num2str(len) ' away from AFp8']);
    
    [AFp4,~] = curveWalk(curve_pts_AFp7AFpzAFp8, 6*len_AFp7AFpzAFp8/8);
    AFp4 = bringPtsToSurf(surf,AFp4);
    %display(['AFp4 = [' num2str(AFp4(1,1)) ' ' num2str(AFp4(1,2)) ' ' num2str(AFp4(1,3)) '], ' num2str(len) ' away from AFp8']);
    
    [AFp6,~] = curveWalk(curve_pts_AFp7AFpzAFp8, 7*len_AFp7AFpzAFp8/8);
    AFp6 = bringPtsToSurf(surf,AFp6);
    %display(['AFp6 = [' num2str(AFp6(1,1)) ' ' num2str(AFp6(1,2)) ' ' num2str(AFp6(1,3)) '], ' num2str(len) ' away from AFp8']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find point along POO7-POOz-POO8 curve 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    plane = plane_equation(POO7, POOz, POO8);
    [curve_pts_POO7POOzPOO8,len_POO7POOzPOO8] = curveGen(POO7, POO8, POOz, plane, surf,dt);
    %display(['POO7-POOz-POO8 curve length: ' num2str(len_POO7POOzPOO8)]);

    [POO5,~] = curveWalk(curve_pts_POO7POOzPOO8, 1*len_POO7POOzPOO8/8);
    POO5 = bringPtsToSurf(surf,POO5);
    %display(['POO5 = [' num2str(POO5(1,1)) ' ' num2str(POO5(1,2)) ' ' num2str(POO5(1,3)) '], ' num2str(len) ' away from POO7']);
        
    [POO3,~] = curveWalk(curve_pts_POO7POOzPOO8, 2*len_POO7POOzPOO8/8);
    POO3 = bringPtsToSurf(surf,POO3);
    %display(['POO3 = [' num2str(POO3(1,1)) ' ' num2str(POO3(1,2)) ' ' num2str(POO3(1,3)) '], ' num2str(len) ' away from POO7']);
 
    [POO1,~] = curveWalk(curve_pts_POO7POOzPOO8, 3*len_POO7POOzPOO8/8);
    POO1 = bringPtsToSurf(surf,POO1);
    %display(['POO1 = [' num2str(POO1(1,1)) ' ' num2str(POO1(1,2)) ' ' num2str(POO1(1,3)) '], ' num2str(len) ' away from POO7']);
       
    [POO2,~] = curveWalk(curve_pts_POO7POOzPOO8, 5*len_POO7POOzPOO8/8);
    POO2 = bringPtsToSurf(surf,POO2);
    %display(['POO2 = [' num2str(POO2(1,1)) ' ' num2str(POO2(1,2)) ' ' num2str(POO2(1,3)) '], ' num2str(len) ' away from POO8']);
      
    [POO4,~] = curveWalk(curve_pts_POO7POOzPOO8, 6*len_POO7POOzPOO8/8);
    POO4 = bringPtsToSurf(surf,POO4);
    %display(['POO4 = [' num2str(POO4(1,1)) ' ' num2str(POO4(1,2)) ' ' num2str(POO4(1,3)) '], ' num2str(len) ' away from POO8']);
    
    [POO6,~] = curveWalk(curve_pts_POO7POOzPOO8, 7*len_POO7POOzPOO8/8);
    POO6 = bringPtsToSurf(surf,POO6);
    %display(['POO6 = [' num2str(POO6(1,1)) ' ' num2str(POO6(1,2)) ' ' num2str(POO6(1,3)) '], ' num2str(len) ' away from POO8']);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find point along NFpz-T9h-OIz curve 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    plane = plane_equation(NFpz, T9h, OIz);
    [curve_pts_NFpzT9hOIz,len_NFpzT9hOIz] = curveGen(NFpz, OIz, T9h, plane, surf, dt);
    %display(['NFpz-T9h-OIz curve length: ' num2str(len_NFpzT9hOIz)]);

    [NFp1h,~] = curveWalk(curve_pts_NFpzT9hOIz, 0.5*len_NFpzT9hOIz/10);
    NFp1h = bringPtsToSurf(surf,NFp1h);
    %display(['NFp1h = [' num2str(NFp1h(1,1)) ' ' num2str(NFp1h(1,2)) ' ' num2str(NFp1h(1,3)) '], ' num2str(len) ' away from NFpz']);
    
    [NFp1,~] = curveWalk(curve_pts_NFpzT9hOIz, 1*len_NFpzT9hOIz/10);
    NFp1 = bringPtsToSurf(surf,NFp1);
    %display(['NFp1 = [' num2str(NFp1(1,1)) ' ' num2str(NFp1(1,2)) ' ' num2str(NFp1(1,3)) '], ' num2str(len) ' away from NFpz']);

    [AFp9h,~] = curveWalk(curve_pts_NFpzT9hOIz, 1.5*len_NFpzT9hOIz/10);
    AFp9h = bringPtsToSurf(surf,AFp9h);
    %display(['AFp9h = [' num2str(AFp9h(1,1)) ' ' num2str(AFp9h(1,2)) ' ' num2str(AFp9h(1,3)) '], ' num2str(len) ' away from NFpz']);
    
    [AF9h,~] = curveWalk(curve_pts_NFpzT9hOIz, 2*len_NFpzT9hOIz/10);
    AF9h = bringPtsToSurf(surf,AF9h);
    %display(['AF9h = [' num2str(AF9h(1,1)) ' ' num2str(AF9h(1,2)) ' ' num2str(AF9h(1,3)) '], ' num2str(len) ' away from NFpz']);
    
    [AFF9h,~] = curveWalk(curve_pts_NFpzT9hOIz, 2.5*len_NFpzT9hOIz/10);
    AFF9h = bringPtsToSurf(surf,AFF9h);
    %display(['AFF9h = [' num2str(AFF9h(1,1)) ' ' num2str(AFF9h(1,2)) ' ' num2str(AFF9h(1,3)) '], ' num2str(len) ' away from NFpz']);
    
    [F9h,~] = curveWalk(curve_pts_NFpzT9hOIz, 3*len_NFpzT9hOIz/10);
    F9h = bringPtsToSurf(surf,F9h);
    %display(['F9h = [' num2str(F9h(1,1)) ' ' num2str(F9h(1,2)) ' ' num2str(F9h(1,3)) '], ' num2str(len) ' away from NFpz']);
    
    [FFT9h,~] = curveWalk(curve_pts_NFpzT9hOIz, 3.5*len_NFpzT9hOIz/10);
    FFT9h = bringPtsToSurf(surf,FFT9h);
    %display(['FFT9h = [' num2str(FFT9h(1,1)) ' ' num2str(FFT9h(1,2)) ' ' num2str(FFT9h(1,3)) '], ' num2str(len) ' away from NFpz']);

    [FT9h,~] = curveWalk(curve_pts_NFpzT9hOIz, 4*len_NFpzT9hOIz/10);
    FT9h = bringPtsToSurf(surf,FT9h);
    %display(['FT9h = [' num2str(FT9h(1,1)) ' ' num2str(FT9h(1,2)) ' ' num2str(FT9h(1,3)) '], ' num2str(len) ' away from NFpz']);
    
    [FTT9h,~] = curveWalk(curve_pts_NFpzT9hOIz, 4.5*len_NFpzT9hOIz/10);
    FTT9h = bringPtsToSurf(surf,FTT9h);
    %display(['FTT9h = [' num2str(FTT9h(1,1)) ' ' num2str(FTT9h(1,2)) ' ' num2str(FTT9h(1,3)) '], ' num2str(len) ' away from NFpz']);
    
    [TTP9h,~] = curveWalk(curve_pts_NFpzT9hOIz, 5.5*len_NFpzT9hOIz/10);
    TTP9h = bringPtsToSurf(surf,TTP9h);
    %display(['TTP9h = [' num2str(TTP9h(1,1)) ' ' num2str(TTP9h(1,2)) ' ' num2str(TTP9h(1,3)) '], ' num2str(len) ' away from NFpz']);
    
    [TP9h,~] = curveWalk(curve_pts_NFpzT9hOIz, 6*len_NFpzT9hOIz/10);
    TP9h = bringPtsToSurf(surf,TP9h);
    %display(['TP9h = [' num2str(TP9h(1,1)) ' ' num2str(TP9h(1,2)) ' ' num2str(TP9h(1,3)) '], ' num2str(len) ' away from NFpz']);
    
    [TPP9h,~] = curveWalk(curve_pts_NFpzT9hOIz, 6.5*len_NFpzT9hOIz/10);
    TPP9h = bringPtsToSurf(surf,TPP9h);
    %display(['TPP9h = [' num2str(TPP9h(1,1)) ' ' num2str(TPP9h(1,2)) ' ' num2str(TPP9h(1,3)) '], ' num2str(len) ' away from NFpz']);
    
    [P9h,~] = curveWalk(curve_pts_NFpzT9hOIz, 7*len_NFpzT9hOIz/10);
    P9h = bringPtsToSurf(surf,P9h);
    %display(['P9h = [' num2str(P9h(1,1)) ' ' num2str(P9h(1,2)) ' ' num2str(P9h(1,3)) '], ' num2str(len) ' away from NFpz']);
    
    [PPO9h,~] = curveWalk(curve_pts_NFpzT9hOIz, 7.5*len_NFpzT9hOIz/10);
    PPO9h = bringPtsToSurf(surf,PPO9h);
    %display(['PPO9h = [' num2str(PPO9h(1,1)) ' ' num2str(PPO9h(1,2)) ' ' num2str(PPO9h(1,3)) '], ' num2str(len) ' away from NFpz']);
    
    [PO9h,~] = curveWalk(curve_pts_NFpzT9hOIz, 8*len_NFpzT9hOIz/10);
    PO9h = bringPtsToSurf(surf,PO9h);
    %display(['PO9h = [' num2str(PO9h(1,1)) ' ' num2str(PO9h(1,2)) ' ' num2str(PO9h(1,3)) '], ' num2str(len) ' away from NFpz']);
    
    [POO9h,~] = curveWalk(curve_pts_NFpzT9hOIz, 8.5*len_NFpzT9hOIz/10);
    POO9h = bringPtsToSurf(surf,POO9h);
    %display(['POO9h = [' num2str(POO9h(1,1)) ' ' num2str(POO9h(1,2)) ' ' num2str(POO9h(1,3)) '], ' num2str(len) ' away from NFpz']);

    [OI1,~] = curveWalk(curve_pts_NFpzT9hOIz, 9*len_NFpzT9hOIz/10);
    OI1 = bringPtsToSurf(surf,OI1);
    %display(['OI1 = [' num2str(OI1(1,1)) ' ' num2str(OI1(1,2)) ' ' num2str(OI1(1,3)) '], ' num2str(len) ' away from NFpz']);
    
    [OI1h,~] = curveWalk(curve_pts_NFpzT9hOIz, 9.5*len_NFpzT9hOIz/10);
    OI1h = bringPtsToSurf(surf,OI1h);
    %display(['OI1h = [' num2str(OI1h(1,1)) ' ' num2str(OI1h(1,2)) ' ' num2str(OI1h(1,3)) '], ' num2str(len) ' away from NFpz']);

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find point along NFpz-T10h-OIz curve 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    plane = plane_equation(NFpz, T10h, OIz);
    [curve_pts_NFpzT10hOIz,len_NFpzT10hOIz] = curveGen(NFpz, OIz, T10h, plane, surf, dt);
    %display(['NFpz-T10h-OIz curve length: ' num2str(len_NFpzT10hOIz)]);

    [NFp2h,~] = curveWalk(curve_pts_NFpzT10hOIz, 0.5*len_NFpzT10hOIz/10);
    NFp2h = bringPtsToSurf(surf,NFp2h);
    %display(['NFp2h = [' num2str(NFp2h(1,1)) ' ' num2str(NFp2h(1,2)) ' ' num2str(NFp2h(1,3)) '], ' num2str(len) ' away from NFpz']);
    
    [NFp2,~] = curveWalk(curve_pts_NFpzT10hOIz, 1*len_NFpzT10hOIz/10);
    NFp2 = bringPtsToSurf(surf,NFp2);
    %display(['NFp2 = [' num2str(NFp2(1,1)) ' ' num2str(NFp2(1,2)) ' ' num2str(NFp2(1,3)) '], ' num2str(len) ' away from NFpz']);

    [AFp10h,~] = curveWalk(curve_pts_NFpzT10hOIz, 1.5*len_NFpzT10hOIz/10);
    AFp10h = bringPtsToSurf(surf,AFp10h);
    %display(['AFp10h = [' num2str(AFp10h(1,1)) ' ' num2str(AFp10h(1,2)) ' ' num2str(AFp10h(1,3)) '], ' num2str(len) ' away from NFpz']);
    
    [AF10h,~] = curveWalk(curve_pts_NFpzT10hOIz, 2*len_NFpzT10hOIz/10);
    AF10h = bringPtsToSurf(surf,AF10h);
    %display(['AF10h = [' num2str(AF10h(1,1)) ' ' num2str(AF10h(1,2)) ' ' num2str(AF10h(1,3)) '], ' num2str(len) ' away from NFpz']);
    
    [AFF10h,~] = curveWalk(curve_pts_NFpzT10hOIz, 2.5*len_NFpzT10hOIz/10);
    AFF10h = bringPtsToSurf(surf,AFF10h);
    %display(['AFF10h = [' num2str(AFF10h(1,1)) ' ' num2str(AFF10h(1,2)) ' ' num2str(AFF10h(1,3)) '], ' num2str(len) ' away from NFpz']);
    
    [F10h,~] = curveWalk(curve_pts_NFpzT10hOIz, 3*len_NFpzT10hOIz/10);
    F10h = bringPtsToSurf(surf,F10h);
    %display(['F10h = [' num2str(F10h(1,1)) ' ' num2str(F10h(1,2)) ' ' num2str(F10h(1,3)) '], ' num2str(len) ' away from NFpz']);
    
    [FFT10h,~] = curveWalk(curve_pts_NFpzT10hOIz, 3.5*len_NFpzT10hOIz/10);
    FFT10h = bringPtsToSurf(surf,FFT10h);
    %display(['FFT10h = [' num2str(FFT10h(1,1)) ' ' num2str(FFT10h(1,2)) ' ' num2str(FFT10h(1,3)) '], ' num2str(len) ' away from NFpz']);

    [FT10h,~] = curveWalk(curve_pts_NFpzT10hOIz, 4*len_NFpzT10hOIz/10);
    FT10h = bringPtsToSurf(surf,FT10h);
    %display(['FT10h = [' num2str(FT10h(1,1)) ' ' num2str(FT10h(1,2)) ' ' num2str(FT10h(1,3)) '], ' num2str(len) ' away from NFpz']);
    
    [FTT10h,~] = curveWalk(curve_pts_NFpzT10hOIz, 4.5*len_NFpzT10hOIz/10);
    FTT10h = bringPtsToSurf(surf,FTT10h);
    %display(['FTT10h = [' num2str(FTT10h(1,1)) ' ' num2str(FTT10h(1,2)) ' ' num2str(FTT10h(1,3)) '], ' num2str(len) ' away from NFpz']);
    
    [TTP10h,~] = curveWalk(curve_pts_NFpzT10hOIz, 5.5*len_NFpzT10hOIz/10);
    TTP10h = bringPtsToSurf(surf,TTP10h);
    %display(['TTP10h = [' num2str(TTP10h(1,1)) ' ' num2str(TTP10h(1,2)) ' ' num2str(TTP10h(1,3)) '], ' num2str(len) ' away from NFpz']);
    
    [TP10h,~] = curveWalk(curve_pts_NFpzT10hOIz, 6*len_NFpzT10hOIz/10);
    TP10h = bringPtsToSurf(surf,TP10h);
    %display(['TP10h = [' num2str(TP10h(1,1)) ' ' num2str(TP10h(1,2)) ' ' num2str(TP10h(1,3)) '], ' num2str(len) ' away from NFpz']);
    
    [TPP10h,~] = curveWalk(curve_pts_NFpzT10hOIz, 6.5*len_NFpzT10hOIz/10);
    TPP10h = bringPtsToSurf(surf,TPP10h);
    %display(['TPP10h = [' num2str(TPP10h(1,1)) ' ' num2str(TPP10h(1,2)) ' ' num2str(TPP10h(1,3)) '], ' num2str(len) ' away from NFpz']);
    
    [P10h,~] = curveWalk(curve_pts_NFpzT10hOIz, 7*len_NFpzT10hOIz/10);
    P10h = bringPtsToSurf(surf,P10h);
    %display(['P10h = [' num2str(P10h(1,1)) ' ' num2str(P10h(1,2)) ' ' num2str(P10h(1,3)) '], ' num2str(len) ' away from NFpz']);
    
    [PPO10h,~] = curveWalk(curve_pts_NFpzT10hOIz, 7.5*len_NFpzT10hOIz/10);
    PPO10h = bringPtsToSurf(surf,PPO10h);
    %display(['PPO10h = [' num2str(PPO10h(1,1)) ' ' num2str(PPO10h(1,2)) ' ' num2str(PPO10h(1,3)) '], ' num2str(len) ' away from NFpz']);
    
    [PO10h,~] = curveWalk(curve_pts_NFpzT10hOIz, 8*len_NFpzT10hOIz/10);
    PO10h = bringPtsToSurf(surf,PO10h);
    %display(['PO10h = [' num2str(PO10h(1,1)) ' ' num2str(PO10h(1,2)) ' ' num2str(PO10h(1,3)) '], ' num2str(len) ' away from NFpz']);
    
    [POO10h,~] = curveWalk(curve_pts_NFpzT10hOIz, 8.5*len_NFpzT10hOIz/10);
    POO10h = bringPtsToSurf(surf,POO10h);
    %display(['POO10h = [' num2str(POO10h(1,1)) ' ' num2str(POO10h(1,2)) ' ' num2str(POO10h(1,3)) '], ' num2str(len) ' away from NFpz']);

    [OI2,~] = curveWalk(curve_pts_NFpzT10hOIz, 9*len_NFpzT10hOIz/10);
    OI2 = bringPtsToSurf(surf,OI2);
    %display(['OI2 = [' num2str(OI2(1,1)) ' ' num2str(OI2(1,2)) ' ' num2str(OI2(1,3)) '], ' num2str(len) ' away from NFpz']);
    
    [OI2h,~] = curveWalk(curve_pts_NFpzT10hOIz, 9.5*len_NFpzT10hOIz/10);
    OI2h = bringPtsToSurf(surf,OI2h);
    %%display(['OI2h = [' num2str(OI2h(1,1)) ' ' num2str(OI2h(1,2)) ' ' num2str(OI2h(1,3)) '], ' num2str(len) ' away from NFpz']);

    % creating the final matrixes
    refpts_10_5 = [...
                    Cz; NFpz; Fpz; AFpz; AFz; AFFz; Fz; FFCz; FCz; FCCz; CCPz; CPz; CPPz; Pz; PPOz; POz; POOz; Oz; OIz; ...
                    T9h; T7; T7h; C5; C5h; C3; C3h; C1; C1h; C2h; C2; C4h; C4; C6h; C6; T8h; T8; T10h; ...
                    Fp1h; Fp1; AFp7; AF7; AFF7; F7; FFT7; FT7; FTT7; TTP7; TP7; TPP7; P7; PPO7; PO7; POO7; O1; O1h; ...
                    Fp2h; Fp2; AFp8; AF8; AFF8; F8; FFT8; FT8; FTT8; TTP8; TP8; TPP8; P8; PPO8; PO8; POO8; O2; O2h; ...
                    F7h; F5; F5h; F3; F3h; F1; F1h; F2h; F2; F4h; F4; F6h; F6; F8h; ...
                    P7h; P5; P5h; P3; P3h; P1; P1h; P2h; P2; P4h; P4; P6h; P6; P8h; ...
                    AF7h; AF5; AF5h; AF3; AF3h; AF1; AF1h; AF2h; AF2; AF4h; AF4; AF6h; AF6; AF8h; ...
                    PO7h; PO5; PO5h; PO3; PO3h; PO1; PO1h; PO2h; PO2; PO4h; PO4; PO6h; PO6; PO8h; ...
                    FT7h; FC5; FC5h; FC3; FC3h; FC1; FC1h; FC2h; FC2; FC4h; FC4; FC6h; FC6; FT8h; ...
                    TP7h; CP5; CP5h; CP3; CP3h; CP1; CP1h; CP2h; CP2; CP4h; CP4; CP6h; CP6; TP8h; ...
                    FTT7h; FCC5; FCC5h; FCC3; FCC3h; FCC1; FCC1h; FCC2h; FCC2; FCC4h; FCC4; FCC6h; FCC6; FTT8h; ...
                    TTP7h; CCP5; CCP5h; CCP3; CCP3h; CCP1; CCP1h; CCP2h; CCP2; CCP4h; CCP4; CCP6h; CCP6; TTP8h; ...
                    FFT7h; FFC5; FFC5h; FFC3; FFC3h; FFC1; FFC1h; FFC2h; FFC2; FFC4h; FFC4; FFC6h; FFC6; FFT8h; ...
                    TPP7h; CPP5; CPP5h; CPP3; CPP3h; CPP1; CPP1h; CPP2h; CPP2; CPP4h; CPP4; CPP6h; CPP6; TPP8h; ...
                    AFF7h; AFF5; AFF5h; AFF3; AFF3h; AFF1; AFF1h; AFF2h; AFF2; AFF4h; AFF4; AFF6h; AFF6; AFF8h; ...
                    PPO7h; PPO5; PPO5h; PPO3; PPO3h; PPO1; PPO1h; PPO2h; PPO2; PPO4h; PPO4; PPO6h; PPO6; PPO8h; ...
                    AFp5; AFp3; AFp1; AFp2; AFp4; AFp6; ...
                    POO5; POO3; POO1; POO2; POO4; POO6; ...
                    NFp1h; NFp1; AFp9h; AF9h; AFF9h; F9h; FFT9h; FT9h; FTT9h; TTP9h; TP9h; TPP9h; P9h; PPO9h; PO9h; POO9h; OI1; OI1h; ...
                    NFp2h; NFp2; AFp10h; AF10h; AFF10h; F10h; FFT10h; FT10h; FTT10h; TTP10h; TP10h; TPP10h; P10h; PPO10h; PO10h; POO10h; OI2; OI2h; ...
                   ];

    refpts_10_5_txt = {...
                        'Cz'; 'NFpz'; 'Fpz'; 'AFpz'; 'AFz'; 'AFFz'; 'Fz'; 'FFCz'; 'FCz'; 'FCCz'; 'CCPz'; 'CPz'; 'CPPz'; 'Pz'; 'PPOz'; 'POz'; 'POOz'; 'Oz'; 'OIz'; ...
                        'T9h'; 'T7'; 'T7h'; 'C5'; 'C5h'; 'C3'; 'C3h'; 'C1'; 'C1h'; 'C2h'; 'C2'; 'C4h'; 'C4'; 'C6h'; 'C6'; 'T8h'; 'T8'; 'T10h'; ...
                        'Fp1h'; 'Fp1'; 'AFp7'; 'AF7'; 'AFF7'; 'F7'; 'FFT7'; 'FT7'; 'FTT7'; 'TTP7'; 'TP7'; 'TPP7'; 'P7'; 'PPO7'; 'PO7'; 'POO7'; 'O1'; 'O1h'; ...
                        'Fp2h'; 'Fp2'; 'AFp8'; 'AF8'; 'AFF8'; 'F8'; 'FFT8'; 'FT8'; 'FTT8'; 'TTP8'; 'TP8'; 'TPP8'; 'P8'; 'PPO8'; 'PO8'; 'POO8'; 'O2'; 'O2h'; ...
                        'F7h'; 'F5'; 'F5h'; 'F3'; 'F3h'; 'F1'; 'F1h'; 'F2h'; 'F2'; 'F4h'; 'F4'; 'F6h'; 'F6'; 'F8h'; ...
                        'P7h'; 'P5'; 'P5h'; 'P3'; 'P3h'; 'P1'; 'P1h'; 'P2h'; 'P2'; 'P4h'; 'P4'; 'P6h'; 'P6'; 'P8h'; ...
                        'AF7h'; 'AF5'; 'AF5h'; 'AF3'; 'AF3h'; 'AF1'; 'AF1h'; 'AF2h'; 'AF2'; 'AF4h'; 'AF4'; 'AF6h'; 'AF6'; 'AF8h'; ...
                        'PO7h'; 'PO5'; 'PO5h'; 'PO3'; 'PO3h'; 'PO1'; 'PO1h'; 'PO2h'; 'PO2'; 'PO4h'; 'PO4'; 'PO6h'; 'PO6'; 'PO8h'; ...
                        'FT7h'; 'FC5'; 'FC5h'; 'FC3'; 'FC3h'; 'FC1'; 'FC1h'; 'FC2h'; 'FC2'; 'FC4h'; 'FC4'; 'FC6h'; 'FC6'; 'FT8h'; ...
                        'TP7h'; 'CP5'; 'CP5h'; 'CP3'; 'CP3h'; 'CP1'; 'CP1h'; 'CP2h'; 'CP2'; 'CP4h'; 'CP4'; 'CP6h'; 'CP6'; 'TP8h'; ...
                        'FTT7h'; 'FCC5'; 'FCC5h'; 'FCC3'; 'FCC3h'; 'FCC1'; 'FCC1h'; 'FCC2h'; 'FCC2'; 'FCC4h'; 'FCC4'; 'FCC6h'; 'FCC6'; 'FTT8h'; ...
                        'TTP7h'; 'CCP5'; 'CCP5h'; 'CCP3'; 'CCP3h'; 'CCP1'; 'CCP1h'; 'CCP2h'; 'CCP2'; 'CCP4h'; 'CCP4'; 'CCP6h'; 'CCP6'; 'TTP8h'; ...
                        'FFT7h'; 'FFC5'; 'FFC5h'; 'FFC3'; 'FFC3h'; 'FFC1'; 'FFC1h'; 'FFC2h'; 'FFC2'; 'FFC4h'; 'FFC4'; 'FFC6h'; 'FFC6'; 'FFT8h'; ...
                        'TPP7h'; 'CPP5'; 'CPP5h'; 'CPP3'; 'CPP3h'; 'CPP1'; 'CPP1h'; 'CPP2h'; 'CPP2'; 'CPP4h'; 'CPP4'; 'CPP6h'; 'CPP6'; 'TPP8h'; ...
                        'AFF7h'; 'AFF5'; 'AFF5h'; 'AFF3'; 'AFF3h'; 'AFF1'; 'AFF1h'; 'AFF2h'; 'AFF2'; 'AFF4h'; 'AFF4'; 'AFF6h'; 'AFF6'; 'AFF8h'; ...
                        'PPO7h'; 'PPO5'; 'PPO5h'; 'PPO3'; 'PPO3h'; 'PPO1'; 'PPO1h'; 'PPO2h'; 'PPO2'; 'PPO4h'; 'PPO4'; 'PPO6h'; 'PPO6'; 'PPO8h'; ...
                        'AFp5'; 'AFp3'; 'AFp1'; 'AFp2'; 'AFp4'; 'AFp6'; ...
                        'POO5'; 'POO3'; 'POO1'; 'POO2'; 'POO4'; 'POO6'; ...
                        'NFp1h'; 'NFp1'; 'AFp9h'; 'AF9h'; 'AFF9h'; 'F9h'; 'FFT9h'; 'FT9h'; 'FTT9h'; 'TTP9h'; 'TP9h'; 'TPP9h'; 'P9h'; 'PPO9h'; 'PO9h'; 'POO9h'; 'OI1'; 'OI1h'; ...
                        'NFp2h'; 'NFp2'; 'AFp10h'; 'AF10h'; 'AFF10h'; 'F10h'; 'FFT10h'; 'FT10h'; 'FTT10h'; 'TTP10h'; 'TP10h'; 'TPP10h'; 'P10h'; 'PPO10h'; 'PO10h'; 'POO10h'; 'OI2'; 'OI2h'; ...
                       };
                   
end

function [curve_seg_best,len_best,max_gap]=curveGen(p1,p2,p3,plane,surf,dt)

%
% USAGE:
%
%    [curve_seg len] = curve_gen(p1, p2, p3, plane, surf, dt)
%
% Input: p1         point from where the curve should start
%        p2         point where the curve should end
%        p3         middle point of the curve
%        plane      4 coefficients of the plane passing through p1, p2 and
%                   p3
%        surf       3D coordinates of the surface of the head
%        dt         threshold setting the distance between the plane and
%                   the points of the surface that will be considered to form the
%                   curve
%
% Output: curve_seg_best   best curve passing through p1, p2 and p3
%         len_best         length of the above curve
%         max_gap          maximum gap between points of the above curve
%
% DESCRIPTION:
%    
%    Find the set of points in surf which lie approximately on the 
%    line of intersection formed by the aurgumant plane and surf. The 
%    output paramter curve_seg is a subset of this set of points 
%    limited to the point lying on the curve segment [p1,p3,p2].
%
% EXAMPLE:
%
%    [curve_pts_NzIz len_NzIz] = curve_gen(Nz, Iz, Cz, plane, surf, .3);
%
% AUTHOR: Jay Dubb (jdubb@nmr.mgh.harvard.edu)
% DATE:   2/9/2010
%
% Modified by: S. Brigadoi 24/04/2013: instead of chosing the best curve as
%                                      the one with less gap between two consecutive 
%                                      points, the initial curve is
%                                      spline interpolated taking 1 every 4 points 
%                                      (to avoid to interpolate also the zig zag between points)
%                                      and the number of points of the curve is increased
%                                      by 5 times. Hence, the curve is now
%                                      more precise. 
%

    if(~exist('dt','var'))
        dt=.6;
    end

    % Generate curve using distance threshhold (dt) from plane of 
    % intersection to points on the surface.   
    [curve_pts] = plane_surf_intersect(plane,surf,dt);
    [curve_seg,~,~] = find_curve_path(p1,p2,p3,curve_pts);

    %display(sprintf('  Curve Init: size=%0.2f, length=%0.2f, threshold=%0.2f, gap size=%0.2f', ...
    %                size(curve_seg,1), len, dt, gap));
      
    % Spline interpolates the curve increasing the number of points 5 times and taking 1 every 4 points.
    % This is done in order not to interpolate the zig zag done by the
    % nodes of the initial curve. Initial and final point of the
    % original curve are included in the interpolated curve
    if mod(size(curve_seg,1)-1,4) ~= 0 
        curve_seg_best = interparc(5*size(curve_seg,1),[curve_seg(1:4:end,1);curve_seg(end,1)],[curve_seg(1:4:end,2);curve_seg(end,2)],[curve_seg(1:4:end,3);curve_seg(end,3)]);
    else
        curve_seg_best = interparc(5*size(curve_seg,1),curve_seg(1:4:end,1),curve_seg(1:4:end,2),curve_seg(1:4:end,3));
    end
    
    % Compute the final length of the curve and the max gap between points
    len = 0;
    max_gap = 0;
    for i = 1:size(curve_seg_best,1)-1
        d = dist3(curve_seg_best(i,:),curve_seg_best(i+1,:));
        len = len+d;
        if(d > max_gap)
            max_gap = d;
        end
    end
    len_best=len;
 
    %display(sprintf('  Final curve: size=%0.2f, length=%0.2f, threshold=%0.2f, gap size=%0.2f', ...
    %                size(curve_seg_best,1), len_best, dt, max_gap));
                
end


function [pf,len_tot] = curveWalk(curve_seg,len)

% USAGE:
%
%    [pf len_tot] = curve_walk(curve_seg, len)
%
% DESCRIPTION:
%    
%    Traverse the curve curve_seg the length of len and return 
%    the coordinates of the point, pf, on which the traversal ends.
%    Also return the distance actually traversed. 
%
%    Do this by summing all the straight-line distances between the 
%    ordered set of points in curve_seg. The first argument curve_seg 
%    must have the property that the array order of each point 
%    corresponds that point's geometric place in the curve. 
%
% EXAMPLE:
%    
%    [Cz  len] = curve_walk(curve_pts_NzIz, 5*len_NzIz/10);
%    
%
% AUTHOR: Jay Dubb (jdubb@nmr.mgh.harvard.edu)
% DATE:   2/9/2010

    len_tot=0;
    for i=1:size(curve_seg,1)-1
        d=dist3(curve_seg(i,:), curve_seg(i+1,:));
        len_tot=len_tot+d;
        if(len_tot>=len)
            break;
        end
    end
    pf=curve_seg(i,:);
    
end

function pts = bringPtsToSurf(surf,pts)

n = size(pts, 1);
for i=1:n
    pts(i,:) = nearestPoint(surf,pts(i,:));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p2_closest,ip2_closest] = nearestPoint(p2,p1)

% AUTHOR: Jay Dubb (jdubb@nmr.mgh.harvard.edu)

m=size(p1,1);

if(~isempty(p2) && ~isempty(p1))
    p2_closest=zeros(m,3);
    ip2_closest=zeros(m,1);
    dmin=zeros(m,1);
    for k=1:m
        d=sqrt((p2(:,1)-p1(k,1)).^2+(p2(:,2)-p1(k,2)).^2+(p2(:,3)-p1(k,3)).^2);
        [dmin(k),ip2_closest(k)]=min(d);
        p2_closest(k,:)=p2(ip2_closest(k),:);
    end
else
    p2_closest=[];
    ip2_closest = 0;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
