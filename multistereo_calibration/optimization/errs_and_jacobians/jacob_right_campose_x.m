function out1 = jacob_right_campose_x(X,Y,ctx,cty,ctz,cwx,cwy,cwz,fx,ptx,pty,ptz,pwx,pwy,pwz,px,rs1_1,rs1_2,rs1_3,rs3_1,rs3_2,rs3_3,s,ts1,ts3)
%JACOB_RIGHT_CAMPOSE_X
%    OUT1 = JACOB_RIGHT_CAMPOSE_X(X,Y,CTX,CTY,CTZ,CWX,CWY,CWZ,FX,PTX,PTY,PTZ,PWX,PWY,PWZ,PX,RS1_1,RS1_2,RS1_3,RS3_1,RS3_2,RS3_3,S,TS1,TS3)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    28-Nov-2018 09:59:41

t2 = cwx.^2;
t3 = cwy.^2;
t4 = cwz.^2;
t5 = t2+t3+t4;
t6 = sqrt(t5);
t7 = cos(t6);
t8 = t7-1.0;
t9 = 1.0./t5;
t10 = sin(t6);
t11 = 1.0./t5.^(3.0./2.0);
t12 = t2+t3;
t13 = 1.0./t5.^2;
t14 = t2+t4;
t15 = 1.0./sqrt(t5);
t16 = t10.*t15;
t17 = t2.*t7.*t9;
t18 = cwx.*cwy.*cwz.*t10.*t11;
t19 = cwx.*cwy.*cwz.*t8.*t13.*2.0;
t20 = cwy.*t2.*t10.*t11;
t21 = cwy.*t2.*t8.*t13.*2.0;
t22 = cwx.*cwz.*t7.*t9;
t23 = cwz.*t2.*t10.*t11;
t24 = cwz.*t2.*t8.*t13.*2.0;
t25 = cwx.*cwy.*t10.*t11;
t26 = t3+t4;
t27 = cwx.*cwy.*t7.*t9;
t28 = cwx.*t10.*t11.*t12;
t29 = cwx.*t8.*t12.*t13.*2.0;
t33 = cwx.*t8.*t9.*2.0;
t30 = t28+t29-t33;
t31 = t2.*t10.*t11;
t32 = cwx.*cwz.*t10.*t11;
t34 = cwx.*t10.*t11.*t14;
t35 = cwx.*t8.*t13.*t14.*2.0;
t36 = -t16-t17+t18+t19+t31;
t37 = ptz.*t36;
t47 = cwy.*t8.*t9;
t38 = t20+t21-t22+t32-t47;
t39 = pty.*t38;
t46 = cwz.*t8.*t9;
t40 = t23+t24-t25+t27-t46;
t41 = ptz.*t40;
t42 = cwx.*t10.*t11.*t26;
t43 = cwx.*t8.*t13.*t26.*2.0;
t44 = t42+t43;
t178 = ptx.*t44;
t45 = t39+t41-t178;
t48 = t23+t24+t25-t27-t46;
t49 = ptx.*t48;
t50 = t16+t17+t18+t19-t31;
t51 = pty.*t50;
t174 = ptz.*t30;
t52 = t49+t51-t174;
t53 = t20+t21+t22-t32-t47;
t54 = ptx.*t53;
t55 = -t33+t34+t35;
t176 = pty.*t55;
t56 = t37+t54-t176;
t57 = pwx.^2;
t58 = pwy.^2;
t59 = pwz.^2;
t60 = t57+t58+t59;
t61 = sqrt(t60);
t62 = sin(t61);
t63 = 1.0./sqrt(t60);
t64 = cos(t61);
t65 = t64-1.0;
t66 = 1.0./t60;
t67 = pwz.*t62.*t63;
t75 = pwx.*pwy.*t65.*t66;
t68 = t67-t75;
t69 = pwy.*t62.*t63;
t70 = pwx.*pwz.*t65.*t66;
t71 = t69+t70;
t72 = t58+t59;
t73 = t65.*t66.*t72;
t74 = t73+1.0;
t76 = t30.*t71;
t77 = t50.*t68;
t78 = t48.*t74;
t79 = t76+t77+t78;
t80 = t55.*t68;
t81 = t36.*t71;
t164 = t53.*t74;
t82 = t80+t81-t164;
t83 = t44.*t74;
t84 = t40.*t71;
t166 = t38.*t68;
t85 = t83+t84-t166;
t86 = t67+t75;
t87 = pwx.*t62.*t63;
t92 = pwy.*pwz.*t65.*t66;
t88 = t87-t92;
t89 = t57+t59;
t90 = t65.*t66.*t89;
t91 = t90+1.0;
t93 = t48.*t86;
t94 = t30.*t88;
t169 = t50.*t91;
t95 = t93+t94-t169;
t96 = t53.*t86;
t97 = t55.*t91;
t171 = t36.*t88;
t98 = t96+t97-t171;
t99 = t40.*t88;
t100 = t44.*t86;
t101 = t38.*t91;
t102 = t99+t100+t101;
t103 = 1.0./s;
t104 = cwz.*t10.*t15;
t105 = cwx.*cwy.*t8.*t9;
t106 = cwy.*t10.*t15;
t107 = cwx.*t10.*t15;
t108 = cwy.*cwz.*t8.*t9;
t109 = t8.*t9.*t12;
t110 = t109+1.0;
t111 = t107-t108;
t112 = cwx.*cwz.*t8.*t9;
t113 = t106+t112;
t114 = t8.*t9.*t14;
t115 = t114+1.0;
t116 = t107+t108;
t117 = t104-t105;
t118 = t104+t105;
t119 = t8.*t9.*t26;
t120 = t119+1.0;
t121 = t106-t112;
t122 = ptx.*t120;
t123 = ptz.*t121;
t181 = pty.*t118;
t124 = ctx+t122+t123-t181;
t125 = rs3_1.*t124;
t126 = pty.*t115;
t127 = ptx.*t117;
t182 = ptz.*t116;
t128 = cty+t126+t127-t182;
t129 = rs3_2.*t128;
t130 = ptz.*t110;
t131 = pty.*t111;
t183 = ptx.*t113;
t132 = ctz+t130+t131-t183;
t133 = rs3_3.*t132;
t134 = t125+t129+t133+ts3;
t135 = t71.*t110;
t136 = t74.*t113;
t184 = t68.*t111;
t137 = t135+t136-t184;
t138 = rs3_3.*t137;
t139 = t68.*t115;
t140 = t71.*t116;
t141 = t74.*t117;
t142 = t139+t140+t141;
t143 = t71.*t121;
t144 = t68.*t118;
t185 = t74.*t120;
t145 = t143+t144-t185;
t146 = rs3_1.*t145;
t186 = rs3_2.*t142;
t147 = t138+t146-t186;
t148 = X.*t147;
t149 = t88.*t110;
t150 = t86.*t113;
t151 = t91.*t111;
t152 = t149+t150+t151;
t153 = t86.*t117;
t154 = t88.*t116;
t187 = t91.*t115;
t155 = t153+t154-t187;
t156 = rs3_2.*t155;
t157 = t86.*t120;
t158 = t91.*t118;
t188 = t88.*t121;
t159 = t157+t158-t188;
t160 = rs3_1.*t159;
t189 = rs3_3.*t152;
t161 = t156+t160-t189;
t162 = Y.*t161;
t244 = t103.*t134;
t163 = t148+t162-t244;
t165 = rs3_2.*t82;
t167 = rs3_1.*t85;
t168 = t165+t167-rs3_3.*t79;
t170 = rs3_3.*t95;
t172 = rs3_2.*t98;
t173 = t170+t172-rs3_1.*t102;
t175 = rs3_3.*t52;
t177 = rs3_2.*t56;
t179 = rs3_1.*t45;
t180 = t175+t177+t179;
t190 = t3.*t10.*t11;
t191 = cwx.*t3.*t10.*t11;
t192 = cwx.*t3.*t8.*t13.*2.0;
t193 = cwy.*cwz.*t10.*t11;
t194 = cwz.*t3.*t10.*t11;
t195 = cwz.*t3.*t8.*t13.*2.0;
t196 = -t25+t27-t46+t194+t195;
t197 = pty.*t196;
t198 = cwy.*t10.*t11.*t12;
t199 = cwy.*t8.*t12.*t13.*2.0;
t203 = cwy.*t8.*t9.*2.0;
t200 = t198+t199-t203;
t201 = t3.*t7.*t9;
t202 = cwy.*cwz.*t7.*t9;
t204 = cwy.*t10.*t11.*t26;
t205 = cwy.*t8.*t13.*t26.*2.0;
t206 = t16+t18+t19-t190+t201;
t207 = ptz.*t206;
t216 = cwx.*t8.*t9;
t208 = t191+t192-t193+t202-t216;
t209 = ptx.*t208;
t210 = t25-t27-t46+t194+t195;
t211 = ptz.*t210;
t212 = cwy.*t10.*t11.*t14;
t213 = cwy.*t8.*t13.*t14.*2.0;
t214 = t212+t213;
t261 = pty.*t214;
t215 = t209+t211-t261;
t217 = -t16+t18+t19+t190-t201;
t218 = ptx.*t217;
t257 = ptz.*t200;
t219 = t197+t218-t257;
t220 = t191+t192+t193-t202-t216;
t221 = pty.*t220;
t222 = -t203+t204+t205;
t259 = ptx.*t222;
t223 = t207+t221-t259;
t224 = t74.*t217;
t225 = t68.*t196;
t226 = t71.*t200;
t227 = t224+t225+t226;
t228 = t71.*t206;
t229 = t74.*t222;
t247 = t68.*t220;
t230 = t228+t229-t247;
t231 = t71.*t210;
t232 = t68.*t214;
t249 = t74.*t208;
t233 = t231+t232-t249;
t234 = t88.*t200;
t235 = t86.*t217;
t252 = t91.*t196;
t236 = t234+t235-t252;
t237 = t86.*t222;
t238 = t88.*t206;
t239 = t91.*t220;
t240 = t237+t238+t239;
t241 = t91.*t214;
t242 = t86.*t208;
t254 = t88.*t210;
t243 = t241+t242-t254;
t245 = 1.0./t163;
t246 = 1.0./t163.^2;
t248 = rs3_1.*t230;
t250 = rs3_2.*t233;
t251 = t248+t250-rs3_3.*t227;
t253 = rs3_3.*t236;
t255 = rs3_2.*t243;
t256 = t253+t255-rs3_1.*t240;
t258 = rs3_3.*t219;
t260 = rs3_1.*t223;
t262 = rs3_2.*t215;
t263 = t258+t260+t262;
t264 = rs1_1.*t124;
t265 = rs1_2.*t128;
t266 = rs1_3.*t132;
t267 = t264+t265+t266+ts1;
t268 = fx.*t267;
t269 = px.*t134;
t270 = t268+t269;
t271 = rs1_3.*t137;
t272 = rs1_1.*t145;
t355 = rs1_2.*t142;
t273 = t271+t272-t355;
t274 = fx.*t273;
t275 = px.*t147;
t276 = t274+t275;
t277 = X.*t276;
t278 = rs1_2.*t155;
t279 = rs1_1.*t159;
t356 = rs1_3.*t152;
t280 = t278+t279-t356;
t281 = fx.*t280;
t282 = px.*t161;
t283 = t281+t282;
t284 = Y.*t283;
t354 = t103.*t270;
t285 = t277+t284-t354;
t286 = t4.*t7.*t9;
t287 = cwx.*t4.*t10.*t11;
t288 = cwx.*t4.*t8.*t13.*2.0;
t289 = cwy.*t4.*t10.*t11;
t290 = cwy.*t4.*t8.*t13.*2.0;
t291 = -t22+t32-t47+t289+t290;
t292 = ptz.*t291;
t293 = cwz.*t10.*t11.*t14;
t294 = cwz.*t8.*t13.*t14.*2.0;
t299 = cwz.*t8.*t9.*2.0;
t295 = t293+t294-t299;
t296 = t4.*t10.*t11;
t297 = -t193+t202-t216+t287+t288;
t298 = ptz.*t297;
t300 = cwz.*t10.*t11.*t26;
t301 = cwz.*t8.*t13.*t26.*2.0;
t302 = -t16+t18+t19-t286+t296;
t303 = pty.*t302;
t304 = t193-t202-t216+t287+t288;
t305 = ptx.*t304;
t306 = t22-t32-t47+t289+t290;
t307 = pty.*t306;
t308 = cwz.*t10.*t11.*t12;
t309 = cwz.*t8.*t12.*t13.*2.0;
t310 = t308+t309;
t351 = ptz.*t310;
t311 = t305+t307-t351;
t312 = t16+t18+t19+t286-t296;
t313 = ptx.*t312;
t347 = pty.*t295;
t314 = t292+t313-t347;
t315 = -t299+t300+t301;
t349 = ptx.*t315;
t316 = t298+t303-t349;
t317 = t71.*t291;
t318 = t68.*t295;
t337 = t74.*t312;
t319 = t317+t318-t337;
t320 = t71.*t297;
t321 = t74.*t315;
t339 = t68.*t302;
t322 = t320+t321-t339;
t323 = t68.*t306;
t324 = t71.*t310;
t325 = t74.*t304;
t326 = t323+t324+t325;
t327 = t86.*t312;
t328 = t91.*t295;
t342 = t88.*t291;
t329 = t327+t328-t342;
t330 = t91.*t302;
t331 = t88.*t297;
t332 = t86.*t315;
t333 = t330+t331+t332;
t334 = t86.*t304;
t335 = t88.*t310;
t344 = t91.*t306;
t336 = t334+t335-t344;
t338 = rs3_2.*t319;
t340 = rs3_1.*t322;
t341 = t338+t340-rs3_3.*t326;
t343 = rs3_2.*t329;
t345 = rs3_3.*t336;
t346 = t343+t345-rs3_1.*t333;
t348 = rs3_2.*t314;
t350 = rs3_1.*t316;
t352 = rs3_3.*t311;
t353 = t348+t350+t352;
out1 = [-t245.*(X.*(px.*t168+fx.*(-rs1_3.*t79+rs1_2.*t82+rs1_1.*t85))+Y.*(px.*t173+fx.*(rs1_3.*t95+rs1_2.*t98-rs1_1.*t102))-t103.*(px.*t180+fx.*(rs1_1.*t45+rs1_3.*t52+rs1_2.*t56)))+t246.*t285.*(X.*t168+Y.*t173-t103.*t180),-t245.*(X.*(px.*t251+fx.*(-rs1_3.*t227+rs1_1.*t230+rs1_2.*t233))+Y.*(px.*t256+fx.*(rs1_3.*t236-rs1_1.*t240+rs1_2.*t243))-t103.*(px.*t263+fx.*(rs1_2.*t215+rs1_3.*t219+rs1_1.*t223)))+t246.*t285.*(X.*t251+Y.*t256-t103.*t263),-t245.*(X.*(px.*t341+fx.*(rs1_2.*t319+rs1_1.*t322-rs1_3.*t326))+Y.*(px.*t346+fx.*(rs1_2.*t329-rs1_1.*t333+rs1_3.*t336))-t103.*(px.*t353+fx.*(rs1_3.*t311+rs1_2.*t314+rs1_1.*t316)))+t246.*t285.*(X.*t341+Y.*t346-t103.*t353),t103.*t245.*(fx.*rs1_1+px.*rs3_1)-rs3_1.*t103.*t246.*t285,t103.*t245.*(fx.*rs1_2+px.*rs3_2)-rs3_2.*t103.*t246.*t285,t103.*t245.*(fx.*rs1_3+px.*rs3_3)-rs3_3.*t103.*t246.*t285];
