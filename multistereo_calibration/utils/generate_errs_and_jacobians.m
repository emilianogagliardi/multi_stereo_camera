function generate_errs_and_jacobians(path)
    txc = sym('ctx'); % inverse camera pose
    tyc = sym('cty');
    tzc = sym('ctz');
    wxc = sym('cwx');
    wyc = sym('cwy');
    wzc = sym('cwz');
    txp = sym('ptx'); % pattern pose
    typ = sym('pty');
    tzp = sym('ptz');
    wxp = sym('pwx');
    wyp = sym('pwy');
    wzp = sym('pwz');
    fx = sym('fx'); % intrinsic parameters
    fy = sym('fy');
    px = sym('px');
    py = sym('py');
    x = sym('x'); % 2D point
    y = sym('y');
    X = sym('X'); % 3D point
    Y = sym('Y');
    Rs = sym('rs', 3);
    ts = sym('ts', [3, 1]);
    s = sym('s');

    K = [fx, 0, px;
         0, fy, py;
         0, 0, 1];

    % rodrigues
    % inverse left camera pose
    thetac = sqrt(wxc^2+wyc^2+wzc^2);
    omegac = [0 -wzc wyc;
             wzc 0 -wxc;
            -wyc wxc 0;];
    Rc = eye(3) + (sin(thetac)/thetac)*omegac + ((1-cos(thetac))/thetac^2)*(omegac*omegac);
    tc = [txc;tyc;tzc];
    
    % pattern pose
    thetap = sqrt(wxp^2+wyp^2+wzp^2);
    omegap = [0 -wzp wyp;
             wzp 0 -wxp;
            -wyp wxp 0;];
    Rp = eye(3) + (sin(thetap)/thetap)*omegap + ((1-cos(thetap))/thetap^2)*(omegap*omegap);
    tp = [txp;typ;tzp];

    % pattern to left camera transformation
    Tl = [Rc, tc; zeros(1, 3), 1] * [Rp, tp; zeros(1, 3), 1];
    Rl = Tl(1:3, 1:3);
    tl = Tl(1:3, 4);
    
    % projection left camera
    projl = K*[Rl(:,1) Rl(:,2) tl]*[X; Y; 1/s];
    projl = projl(1:2, :) ./projl(3, :);
    
    % error and jacobian left camera
    ol = [x; y] - projl;
    jlc = jacobian(ol, [wxc, wyc, wzc, txc, tyc, tzc]);
    jlp = jacobian(ol, [wxp, wyp, wzp, txp, typ, tzp]);
    jls = jacobian(ol, s);
    
    % pattern to right camera transformation
    Tr = [Rs, ts; zeros(1, 3), 1] * Tl;
    Rr = Tr(1:3, 1:3);
    tr = Tr(1:3, 4);
    
    % projection right camera
    projr = K*[Rr(:,1) Rr(:,2) tr]*[X; Y; 1/s];
    projr = projr(1:2, :) ./projr(3, :);
    
    % error and jacobian left camera
    or = [x; y] - projr;
    jrc = jacobian(or, [wxc, wyc, wzc, txc, tyc, tzc]);
    jrp = jacobian(or, [wxp, wyp, wzp, txp, typ, tzp]);
    jrs = jacobian(or, s);
    
    % save the functions
    matlabFunction(ol(1), 'File', [path, '/err_left_x.m']);
    matlabFunction(ol(2), 'File', [path, '/err_left_y.m']);
    matlabFunction(or(1), 'File', [path, '/err_right_x.m']);
    matlabFunction(or(2), 'File', [path, '/err_right_y.m']);
    matlabFunction(jlc(1, :), 'File', [path, '/jacob_left_campose_x.m']);
    matlabFunction(jlc(2, :), 'File', [path, '/jacob_left_campose_y.m']);
    matlabFunction(jrc(1, :), 'File', [path, '/jacob_right_campose_x.m']);
    matlabFunction(jrc(2, :), 'File', [path, '/jacob_right_campose_y.m']);
    matlabFunction(jlp(1, :), 'File', [path, '/jacob_left_patpose_x.m']);
    matlabFunction(jlp(2, :), 'File', [path, '/jacob_left_patpose_y.m']);
    matlabFunction(jrp(1, :), 'File', [path, '/jacob_right_patpose_x.m']);
    matlabFunction(jrp(2, :), 'File', [path, '/jacob_right_patpose_y.m']);
    matlabFunction(jls(1, :), 'File', [path, '/jacob_left_s_x.m']);
    matlabFunction(jls(2, :), 'File', [path, '/jacob_left_s_y.m']);
    matlabFunction(jrs(1, :), 'File', [path, '/jacob_right_s_x.m']);
    matlabFunction(jrs(2, :), 'File', [path, '/jacob_right_s_y.m']);
end
