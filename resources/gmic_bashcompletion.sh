#
#  Bash completion rules for 'gmic'.
#
# This file has been generated automatically.
# Do not edit!
#
# This file should be copied/renamed in '/usr/share/bash-completion/completions/gmic'.
#

_gmic()
{
    local cur prev opts coms
    if type -t _init_completion >/dev/null; then
        _init_completion -n = || return
    else
        COMPREPLY=()
        cur="${COMP_WORDS[COMP_CWORD]}"
        prev="${COMP_WORDS[COMP_CWORD-1]}"
    fi
    coms="!= % & *3d * +3d + -3d - /3d / << <= < == => = >= >> > ^ a abs ac acos acosh add3d add adjust_colors alert and animate3d animate ap apc apo append append_tiles apply_camera3d apply_camera apply_channels apply_curve apply_files apply_gamma apply_matrix3d apply_parallel apply_parallel_channels apply_parallel_overlap apply_scales apply_tiles apply_timeout apply_video area area_fg arg0 arg2img arg2var arg argmax argmaxabs argmin argminabs array3d array array_fade array_mirror array_random arrow3d arrow asin asinh at at_line at_quadrangle atan2 atan atanh autocrop autocrop_components autocrop_coords autocrop_seq autoindex average_files average_vectors average_video axes3d axes b balance_gamma ball bandpass barycenter base642img base642uint8 basename bayer2rgb bilateral bin2dec bin blend blend blend_edges blend_fade blend_median blend_seamless blur blur_angular blur_bloom blur_linear blur_radial blur_selective blur_x blur_xy blur_xyz blur_y blur_z boundingbox3d box3d boxfilter boxfitting break brushify bsl bsr bump2normal c3d c camera cartoon cast center3d channels check3d check chessboard cie1931 circle3d circle circles3d close_binary closing closing_circ clut cmy2rgb cmyk2rgb col3d color3d color_ellipses colorblind colorcube3d colormap columns command complex2polar compose_channels compose_freq compress_clut compress_rle cone3d continue convolve convolve_fft correlate cos cosh covariance_vectors cracks crop cross_correlation cubes3d cubism cumulate cup3d cursor curvature cut cylinder3d d0 d2d d3d d da db3d dc dct deblur deblur_goldmeinel deblur_richardsonlucy debug dec2bin dec2hex dec2oct dec2str dec decompress_clut decompress_clut_pde decompress_clut_rbf decompress_rle deconvolve_fft deform deg2rad deinterlace delaunay3d delaunay delete deltaE demos denoise denoise_cnn denoise_haar denoise_patchpca deriche detect_skin dfft dg dh diagonal diffusiontensors dijkstra dilate dilate_circ dilate_oct dilate_threshold direction2rgb discard displacement display0 display2d display3d display display_array display_camera display_fft display_graph display_histogram display_parallel0 display_parallel display_parametric display_polar display_quiver display_rgba display_tensors display_warp distance distribution3d ditheredbw div3d div div_complex divergence do dog done double3d dp0 dp dq draw_whirl drawing drgba drop_shadow dt dw e echo echo_file edgels edges eigen2tensor eigen elevate elevation3d elif ellipse ellipsionism else empty3d endian eq equalize equirectangular2nadirzenith erf erode erode_circ erode_oct erode_threshold error euclidean2polar eval exec exec_out exp expand_x expand_xy expand_xyz expand_y expand_z extract extract_region extrude3d eye f3d f fact fade_diamond fade_files fade_linear fade_radial fade_video fade_x fade_y fade_z fc fft fftpolar fi fibonacci file_mv file_rand filename files2img files2video files fill fill_color fire_edges fisheye fitratio_wh fitscreen flood flower focale3d fontchart for foreach fps fractalize frame frame_blur frame_cube frame_fuzzy frame_painting frame_pattern frame_round frame_seamless frame_x frame_xy frame_xyz frame_y function1d g gaussian gaussians3d gcd ge glow gmd2ascii gmd2html gmic3d gradient2rgb gradient gradient_norm gradient_orientation graph grid gt guided gyroid3d h haar halftone hardsketchbw hcy2rgb hearts heat_flow help hessian hex2dec hex2img8 hex2img hex2str hex histogram3d histogram histogram_cumul histogram_nd histogram_pointwise hough houghsketchbw hsi2rgb hsi82rgb hsl2rgb hsl82rgb hsv2rgb hsv82rgb i idct identity iee if ifft ifftpolar ig ihaar ilaplacian image6cube3d image imageblocks3d imagecube3d imagegrid imagegrid_hexagonal imagegrid_triangular imageplane3d imagepyramid3d imagerubik3d imagesphere3d img2ascii img2base64 img2hex img2patches img2str img2text img82hex index inn inpaint inpaint_flow inpaint_holes inpaint_matchpatch inpaint_morpho inpaint_pde input input_565 input_cached input_csv input_cube input_flo input_glob input_gpl input_obj input_text inrange int2rgb invert ipremula ir is_3d is_change is_ext is_half is_image_arg is_macos is_pattern is_percent is_videofilename is_windows isoline3d isophotes isosurface3d it j3d j jzazbz2rgb jzazbz2xyz k kaleidoscope keep keep_named kn kuwahara l3d l laar lab2lch lab2rgb lab2srgb lab2xyz lab82rgb lab82srgb label3d label label_fg label_points3d laplacian lathe3d lch2lab lch2rgb lch82rgb le lic light3d light_patch light_relief lightness lightrays line3d line linearize_tiles linify lissajous3d local log10 log2 log lorem lt luminance lut_contrast m* m/ m3d m mad mandelbrot map map_clut map_sphere map_sprites map_tones map_tones_fast marble matchpatch math_lib max max_d max_h max_patch max_s max_w max_wh max_whd max_whds maxabs maze maze_mask md3d mdiv meancurvature_flow med median median_files median_vectors median_video meigen min min_d min_h min_patch min_s min_w min_wh min_whd min_whds minabs minimal_path mirror mix_channels mix_rgb mmul mod mode3d moded3d montage morph morph_files morph_rbf morph_video mosaic move mproj mse mse_matrix mul3d mul mul_channels mul_complex mutex mv n3d n nadirzenith2equirectangular name named negate neq network newton_fractal nlmeans nlmeans_core nm nmd nn_check_layer nn_init nn_layer_add nn_layer_append nn_layer_avgpool2d nn_layer_batchnorm nn_layer_clone nn_layer_conv2d nn_layer_conv2dbnnl nn_layer_conv2dnl nn_layer_crop nn_layer_fc nn_layer_fcbnnl nn_layer_fcnl nn_layer_input nn_layer_maxpool2d nn_layer_nl nn_layer_rename nn_layer_resconv2d nn_layer_resconv2dnl nn_layer_reshape nn_layer_resize nn_layer_run nn_layer_split nn_lib nn_load nn_loss_bce nn_loss_mse nn_save nn_trainer noarg noise noise_hurl noise_perlin noise_poissondisk norm normalize3d normalize normalize_filename normalize_l2 normalize_local normalize_sum normalized_cross_correlation normp not o3d o object3d oct2dec oct oklab2rgb old_photo on oneminus onfail op opacity3d opening opening_circ or orientation orthogonalize ot otsu output output_565 output_cube output_flo output_ggr output_gmz output_obj output_text outputn outputp outputw outputx ow ox p3d p pack pack_sprites padint palette parallel parametric3d parse_cli parse_gmd parse_gui pass patches2img patches path_cache path_current path_gimp path_tmp pca_patch3d pde_flow pencilbw percentile periodize_poisson permute peronamalik_flow phase_correlation piechart pixelize pixelsort plane3d plasma plot2value plot point3d point pointcloud3d pointcloud polar2complex polar2euclidean polaroid polka_dots polygon polygonize portrait pose3d poster_edges poster_hope pow premula primitives3d print progress projections3d pseudogray psnr psnr_matrix puzzle pyramid3d q quadrangle3d quadratize_tiles quantize quantize_area quit quiver r2din r2dout r2dx r2dy r3d r3din r3dout r3dx r3dy r3dz r rad2deg raindrops rand random3d random_pattern rbf rectangle red_eye register_nonrigid register_rigid remove remove_copymark remove_duplicates remove_empty remove_hotpixels remove_named remove_opacity remove_pixels repeat replace replace_color replace_inf replace_nan replace_naninf replace_seq replace_str reset resize2din resize2dout resize2dx resize2dy resize3din resize3dout resize3dx resize3dy resize3dz resize resize_as_image resize_mn resize_pow2 resize_ratio2d retinex return reverse3d reverse rgb2bayer rgb2cmy rgb2cmyk rgb2hcy rgb2hsi8 rgb2hsi rgb2hsl8 rgb2hsl rgb2hsv8 rgb2hsv rgb2int rgb2jzazbz rgb2lab8 rgb2lab rgb2lch8 rgb2lch rgb2luv rgb2oklab rgb2ryb rgb2srgb rgb2xyz8 rgb2xyz rgb2ycbcr rgb2yiq8 rgb2yiq rgb2yuv8 rgb2yuv rgb rgba ri ripple rm rmn rodilius rol rolling_guidance ror rorschach rotate3d rotate rotate_tileable rotate_tiles rotation3d rotoidoscope round roundify rows rprogress rr2d run rv3d rv ryb2rgb s3d s sample scale2x scale3x scale_dcci2x scanlines screen seamcarve segment_watershed select select_color sepia serialize set sh shade_stripes shadow_patch shape2bump shape_circle shape_cupid shape_diamond shape_dragon shape_fern shape_gear shape_heart shape_polygon shape_snowflake shape_star shared sharpen shell_cols shift shift_tiles shrink_x shrink_xy shrink_xyz shrink_y shrink_z sierpinski3d sierpinski sign sin sinc sinh size3d size_value skeleton3d skeleton sketchbw skip sl3d slic slices smooth snapshot3d solarize solidify solve solve_poisson sort sort_list sp specl3d specs3d sphere3d spherical3d spherize spiralbw spline3d spline split3d split split_colors split_details split_freq split_opacity split_tiles sponge spread sprite3d sprites3d sqr sqrt srand srgb2lab8 srgb2lab srgb2rgb ss3d ssd_patch ssim ssim_matrix stained_glass star3d stars status std_noise stencil stencilbw store str2hex str strcapitalize strcasevar strcontains streamline3d stripes_y strlen strlowercase strreplace structuretensors struppercase strvar strver stylize sub3d sub sub_alpha superformula3d surfels3d svd symmetrize syntexturize syntexturize_matchpatch t3d t tan tanh taquin tensors3d testimage2d tetraedron_shade tetris text3d text text_outline text_pointcloud3d texturize3d texturize_canvas texturize_paper thickline thinning threshold tic tixy to to_a to_clutname to_color to_colormode to_gray to_graya to_pseudogray to_rgb to_rgba toc tones topographic_map torus3d transfer_histogram transfer_pca transfer_rgb transform_polar transition3d transition transpose triangle3d triangle_shade trisolve truchet tsp tunnel turbulence tv_flow twirl u uint82base64 um uncommand undistort uniform_distribution unroll unserialize unsharp unsharp_octave up update upscale_smart v vanvliet variance_patch vector2tensor verbose version video2files vignette volume3d voronoi w wait warhol warn warp warp_patch warp_perspective warp_rbf water watermark_fourier watermark_visible watershed wave weave weird3d while whirls wind window x x_2048 x_blobs x_bouncing x_color_curves x_colorize x_connect4 x_crop x_cut x_fire x_fireworks x_fisheye x_fourier x_grab_color x_hanoi x_histogram x_hough x_jawbreaker x_landscape x_life x_light x_mandelbrot x_mask_color x_metaballs3d x_minesweeper x_minimal_path x_morph x_pacman x_paint x_plasma x_quantize_rgb x_reflection3d x_rubber3d x_segment x_select_color x_select_function1d x_select_palette x_shadebobs x_spline x_starfield3d x_tetris x_threshold x_tictactoe x_warp x_waves x_whirl xo xor xyz2jzazbz xyz2lab xyz2rgb xyz82rgb xz y ycbcr2rgb yinyang yiq2rgb yiq82rgb yuv2rgb yuv82rgb z zoom | }"
    opts=$(echo "$coms" | sed "s: \([^ ]\+\): \1 -\1 \+\1:g")

    case "${prev}" in
        "!=" | "-!=" | "+!=")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "%" | "-%" | "+%")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "&" | "-&" | "+&")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "*3d" | "-*3d" | "+*3d")
            COMPREPLY=( $(compgen -W "factor factor_x,factor_y,_factor_z") ); return 0;;
        "*" | "-*" | "+*")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "+3d" | "-+3d" | "++3d")
            COMPREPLY=( $(compgen -W "tx,_ty,_tz [object3d] (no_arg)") ); return 0;;
        "+" | "-+" | "++")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "-3d" | "--3d" | "+-3d")
            COMPREPLY=( $(compgen -W "> tx,_ty,_tz") ); return 0;;
        "-" | "--" | "+-")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "/3d" | "-/3d" | "+/3d")
            COMPREPLY=( $(compgen -W "factor factor_x,factor_y,_factor_z") ); return 0;;
        "/" | "-/" | "+/")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "<<" | "-<<" | "+<<")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "<=" | "-<=" | "+<=")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "<" | "-<" | "+<")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "==" | "-==" | "+==")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "=>" | "-=>" | "+=>")
            COMPREPLY=( $(compgen -W "> \"name1\",\"name2\",...") ); return 0;;
        "=" | "-=" | "+=")
            COMPREPLY=( $(compgen -W "> value,_x[%],_y[%],_z[%],_c[%]") ); return 0;;
        ">=" | "->=" | "+>=")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        ">>" | "->>" | "+>>")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        ">" | "->" | "+>")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "^" | "-^" | "+^")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "a" | "-a" | "+a")
            COMPREPLY=( $(compgen -W "[image],axis,_centering axis,_centering") ); return 0;;
        "ac" | "-ac" | "+ac")
            COMPREPLY=( $(compgen -W "> \"command\",color_channels,_value_action={_0=none_|_1=cut_|_2=normalize_}") ); return 0;;
        "add3d" | "-add3d" | "+add3d")
            COMPREPLY=( $(compgen -W "tx,_ty,_tz [object3d] (no_arg)") ); return 0;;
        "add" | "-add" | "+add")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "adjust_colors" | "-adjust_colors" | "+adjust_colors")
            COMPREPLY=( $(compgen -W "> -100<=_brightness<=100,-100<=_contrast<=100,-100<=_gamma<=100,-100<=_hue_shift<=100,-100<=_saturation<=100,_value_min,_value_max") ); return 0;;
        "alert" | "-alert" | "+alert")
            COMPREPLY=( $(compgen -W "> _title,_message,_label_button1,_label_button2,...") ); return 0;;
        "and" | "-and" | "+and")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "animate3d" | "-animate3d" | "+animate3d")
            COMPREPLY=( $(compgen -W "> nb_frames>0,_step_angle_x,_step_angle_y,_step_angle_z,_zoom_factor,0<=_fake_shadow_level<=100,_[background]") ); return 0;;
        "animate" | "-animate" | "+animate")
            COMPREPLY=( $(compgen -W "filter_name,\"param1_start,...,paramN_start\",\"param1_end,...,paramN_end\",nb_frames>=0,_output_frames={_0_|_1_},_output_filename delay>0,_back_and_forth={_0_|_1_}") ); return 0;;
        "ap" | "-ap" | "+ap")
            COMPREPLY=( $(compgen -W "> \"command\"") ); return 0;;
        "apc" | "-apc" | "+apc")
            COMPREPLY=( $(compgen -W "> \"command\"") ); return 0;;
        "apo" | "-apo" | "+apo")
            COMPREPLY=( $(compgen -W "> \"command\",overlap[%],nb_threads={_0=auto_|_1_|_2_|_4_|_8_|_16_}") ); return 0;;
        "append" | "-append" | "+append")
            COMPREPLY=( $(compgen -W "[image],axis,_centering axis,_centering") ); return 0;;
        "append_tiles" | "-append_tiles" | "+append_tiles")
            COMPREPLY=( $(compgen -W "> _M>=0,_N>=0,0<=_centering_x<=1,0<=_centering_y<=1") ); return 0;;
        "apply_camera3d" | "-apply_camera3d" | "+apply_camera3d")
            COMPREPLY=( $(compgen -W "> pos_x,pos_y,pos_z,target_x,target_y,target_z,up_x,up_y,up_z") ); return 0;;
        "apply_camera" | "-apply_camera" | "+apply_camera")
            COMPREPLY=( $(compgen -W "> _\"command\",_camera_index>=0,_skip_frames>=0,_output_filename") ); return 0;;
        "apply_channels" | "-apply_channels" | "+apply_channels")
            COMPREPLY=( $(compgen -W "> \"command\",color_channels,_value_action={_0=none_|_1=cut_|_2=normalize_}") ); return 0;;
        "apply_curve" | "-apply_curve" | "+apply_curve")
            COMPREPLY=( $(compgen -W "> 0<=smoothness<=1,x0,y0,x1,y1,x2,y2,...,xN,yN") ); return 0;;
        "apply_files" | "-apply_files" | "+apply_files")
            COMPREPLY=( $(compgen -W "> \"filename_pattern\",_\"command\",_first_frame>=0,_last_frame={_>=0_|_-1=last_},_frame_step>=1,_output_filename") ); return 0;;
        "apply_gamma" | "-apply_gamma" | "+apply_gamma")
            COMPREPLY=( $(compgen -W "> gamma>=0") ); return 0;;
        "apply_matrix3d" | "-apply_matrix3d" | "+apply_matrix3d")
            COMPREPLY=( $(compgen -W "> a11,a12,a13,...,a31,a32,a33") ); return 0;;
        "apply_parallel" | "-apply_parallel" | "+apply_parallel")
            COMPREPLY=( $(compgen -W "> \"command\"") ); return 0;;
        "apply_parallel_channels" | "-apply_parallel_channels" | "+apply_parallel_channels")
            COMPREPLY=( $(compgen -W "> \"command\"") ); return 0;;
        "apply_parallel_overlap" | "-apply_parallel_overlap" | "+apply_parallel_overlap")
            COMPREPLY=( $(compgen -W "> \"command\",overlap[%],nb_threads={_0=auto_|_1_|_2_|_4_|_8_|_16_}") ); return 0;;
        "apply_scales" | "-apply_scales" | "+apply_scales")
            COMPREPLY=( $(compgen -W "> \"command\",number_of_scales>0,_min_scale[%]>=0,_max_scale[%]>=0,_scale_gamma>0,_interpolation") ); return 0;;
        "apply_tiles" | "-apply_tiles" | "+apply_tiles")
            COMPREPLY=( $(compgen -W "> \"command\",_tile_width[%]>0,_tile_height[%]>0,_tile_depth[%]>0,_overlap_width[%]>=0,_overlap_height[%]>=0,_overlap_depth[%]>=0,_boundary_conditions={_0=dirichlet_|_1=neumann_|_2=periodic_|_3=mirror_}") ); return 0;;
        "apply_timeout" | "-apply_timeout" | "+apply_timeout")
            COMPREPLY=( $(compgen -W "> \"command\",_timeout={_0=no_timeout_|_>0=with_specified_timeout_(in_seconds)_}") ); return 0;;
        "apply_video" | "-apply_video" | "+apply_video")
            COMPREPLY=( $(compgen -W "> video_filename,_\"command\",_first_frame>=0,_last_frame={_>=0_|_-1=last_},_frame_step>=1,_output_filename") ); return 0;;
        "area" | "-area" | "+area")
            COMPREPLY=( $(compgen -W "> tolerance>=0,is_high_connectivity={_0_|_1_}") ); return 0;;
        "area_fg" | "-area_fg" | "+area_fg")
            COMPREPLY=( $(compgen -W "> tolerance>=0,is_high_connectivity={_0_|_1_}") ); return 0;;
        "arg0" | "-arg0" | "+arg0")
            COMPREPLY=( $(compgen -W "> n>=0,_arg0,...,_argN") ); return 0;;
        "arg2img" | "-arg2img" | "+arg2img")
            COMPREPLY=( $(compgen -W "> argument_1,...,argument_N") ); return 0;;
        "arg2var" | "-arg2var" | "+arg2var")
            COMPREPLY=( $(compgen -W "> variable_name,argument_1,...,argument_N") ); return 0;;
        "arg" | "-arg" | "+arg")
            COMPREPLY=( $(compgen -W "> n>=1,_arg1,...,_argN") ); return 0;;
        "array3d" | "-array3d" | "+array3d")
            COMPREPLY=( $(compgen -W "> size_x>=1,_size_y>=1,_size_z>=1,_offset_x[%],_offset_y[%],_offset_y[%]") ); return 0;;
        "array" | "-array" | "+array")
            COMPREPLY=( $(compgen -W "> M>0,_N>0,_expand_type={_0=min_|_1=max_|_2=all_}") ); return 0;;
        "array_fade" | "-array_fade" | "+array_fade")
            COMPREPLY=( $(compgen -W "> M>0,_N>0,0<=_fade_start<=100,0<=_fade_end<=100,_expand_type={0=min_|_1=max_|_2=all}") ); return 0;;
        "array_mirror" | "-array_mirror" | "+array_mirror")
            COMPREPLY=( $(compgen -W "> N>=0,_dir={_0=x_|_1=y_|_2=xy_|_3=tri-xy_},_expand_type={_0_|_1_}") ); return 0;;
        "array_random" | "-array_random" | "+array_random")
            COMPREPLY=( $(compgen -W "> Ms>0,_Ns>0,_Md>0,_Nd>0") ); return 0;;
        "arrow3d" | "-arrow3d" | "+arrow3d")
            COMPREPLY=( $(compgen -W "> x0,y0,z0,x1,y1,z1,_radius[%]>=0,_head_length[%]>=0,_head_radius[%]>=0") ); return 0;;
        "arrow" | "-arrow" | "+arrow")
            COMPREPLY=( $(compgen -W "> x0[%],y0[%],x1[%],y1[%],_thickness[%]>=0,_head_length[%]>=0,_head_thickness[%]>=0,_opacity,_pattern,_color1,...") ); return 0;;
        "at" | "-at" | "+at")
            COMPREPLY=( $(compgen -W "> \"command\",_tile_width[%]>0,_tile_height[%]>0,_tile_depth[%]>0,_overlap_width[%]>=0,_overlap_height[%]>=0,_overlap_depth[%]>=0,_boundary_conditions={_0=dirichlet_|_1=neumann_|_2=periodic_|_3=mirror_}") ); return 0;;
        "at_line" | "-at_line" | "+at_line")
            COMPREPLY=( $(compgen -W "> x0[%],y0[%],z0[%],x1[%],y1[%],z1[%]") ); return 0;;
        "at_quadrangle" | "-at_quadrangle" | "+at_quadrangle")
            COMPREPLY=( $(compgen -W "x0[%],y0[%],x1[%],y1[%],x2[%],y2[%],x3[%],y3[%],_interpolation,_boundary_conditions x0[%],y0[%],z0[%],x1[%],y1[%],z1[%],x2[%],y2[%],z2[%],x3[%],y3[%],z3[%],_interpolation,_boundary_conditions") ); return 0;;
        "atan2" | "-atan2" | "+atan2")
            COMPREPLY=( $(compgen -W "> [x_argument]") ); return 0;;
        "autocrop" | "-autocrop" | "+autocrop")
            COMPREPLY=( $(compgen -W "value1,value2,... (no_arg)") ); return 0;;
        "autocrop_components" | "-autocrop_components" | "+autocrop_components")
            COMPREPLY=( $(compgen -W "> _threshold[%],_min_area[%]>=0,_is_high_connectivity={_0_|_1_},_output_type={_0=crop_|_1=segmentation_|_2=coordinates_}") ); return 0;;
        "autocrop_coords" | "-autocrop_coords" | "+autocrop_coords")
            COMPREPLY=( $(compgen -W "> value1,value2,..._|_auto") ); return 0;;
        "autocrop_seq" | "-autocrop_seq" | "+autocrop_seq")
            COMPREPLY=( $(compgen -W "> value1,value2,..._|_auto") ); return 0;;
        "autoindex" | "-autoindex" | "+autoindex")
            COMPREPLY=( $(compgen -W "> nb_colors>0,0<=_dithering<=1,_method={_0=median-cut_|_1=k-means_}") ); return 0;;
        "average_files" | "-average_files" | "+average_files")
            COMPREPLY=( $(compgen -W "> \"filename_pattern\",_first_frame>=0,_last_frame={_>=0_|_-1=last_},_frame_step>=1,_output_filename") ); return 0;;
        "average_video" | "-average_video" | "+average_video")
            COMPREPLY=( $(compgen -W "> video_filename,_first_frame>=0,_last_frame={_>=0_|_-1=last_},_frame_step>=1,_output_filename") ); return 0;;
        "axes3d" | "-axes3d" | "+axes3d")
            COMPREPLY=( $(compgen -W "> _size_x,_size_y,_size_z,_font_size>0,_label_x,_label_y,_label_z,_is_origin={_0=no_|_1=yes_}") ); return 0;;
        "axes" | "-axes" | "+axes")
            COMPREPLY=( $(compgen -W "> x0,x1,y0,y1,_font_height>=0,_opacity,_pattern,_color1,...") ); return 0;;
        "b" | "-b" | "+b")
            COMPREPLY=( $(compgen -W "std_deviation>=0[%],_boundary_conditions,_kernel axes,std_deviation>=0[%],_boundary_conditions,_kernel") ); return 0;;
        "balance_gamma" | "-balance_gamma" | "+balance_gamma")
            COMPREPLY=( $(compgen -W "> _ref_color1,...") ); return 0;;
        "ball" | "-ball" | "+ball")
            COMPREPLY=( $(compgen -W "> _size>0,__R,_G,_B,0<=_specular_light<=8,0<=_specular_size<=8,_shadow>=0") ); return 0;;
        "bandpass" | "-bandpass" | "+bandpass")
            COMPREPLY=( $(compgen -W "> _min_freq[%],_max_freq[%]") ); return 0;;
        "base642img" | "-base642img" | "+base642img")
            COMPREPLY=( $(compgen -W "> \"base64_string\"") ); return 0;;
        "base642uint8" | "-base642uint8" | "+base642uint8")
            COMPREPLY=( $(compgen -W "> \"base64_string\"") ); return 0;;
        "basename" | "-basename" | "+basename")
            COMPREPLY=( $(compgen -W "> file_path,_variable_name_for_folder") ); return 0;;
        "bayer2rgb" | "-bayer2rgb" | "+bayer2rgb")
            COMPREPLY=( $(compgen -W "> _GM_smoothness,_RB_smoothness1,_RB_smoothness2") ); return 0;;
        "bilateral" | "-bilateral" | "+bilateral")
            COMPREPLY=( $(compgen -W "[guide],std_deviation_s[%]>=0,std_deviation_r[%]>=0,_sampling_s>=0,_sampling_r>=0 std_deviation_s[%]>=0,std_deviation_r[%]>=0,_sampling_s>=0,_sampling_r>=0") ); return 0;;
        "bin2dec" | "-bin2dec" | "+bin2dec")
            COMPREPLY=( $(compgen -W "> binary_int1,...") ); return 0;;
        "bin" | "-bin" | "+bin")
            COMPREPLY=( $(compgen -W "> binary_int1,...") ); return 0;;
        "blend" | "-blend" | "+blend")
            COMPREPLY=( $(compgen -W "[layer],blending_mode,_opacity[%],_selection_is={_0=base-layers_|_1=top-layers_} blending_mode,_opacity[%]") ); return 0;;
        "blend" | "-blend" | "+blend")
            COMPREPLY=( $(compgen -W "[layer],blending_mode,_opacity[%],_selection_is={_0=base-layers_|_1=top-layers_} blending_mode,_opacity[%]") ); return 0;;
        "blend_edges" | "-blend_edges" | "+blend_edges")
            COMPREPLY=( $(compgen -W "> smoothness[%]>=0") ); return 0;;
        "blend_fade" | "-blend_fade" | "+blend_fade")
            COMPREPLY=( $(compgen -W "> [fading_shape]") ); return 0;;
        "blend_seamless" | "-blend_seamless" | "+blend_seamless")
            COMPREPLY=( $(compgen -W "> _is_mixed_mode={_0_|_1_},_inner_fading[%]>=0,_outer_fading[%]>=0") ); return 0;;
        "blur" | "-blur" | "+blur")
            COMPREPLY=( $(compgen -W "std_deviation>=0[%],_boundary_conditions,_kernel axes,std_deviation>=0[%],_boundary_conditions,_kernel") ); return 0;;
        "blur_angular" | "-blur_angular" | "+blur_angular")
            COMPREPLY=( $(compgen -W "> amplitude[%],_center_x[%],_center_y[%]") ); return 0;;
        "blur_bloom" | "-blur_bloom" | "+blur_bloom")
            COMPREPLY=( $(compgen -W "> _amplitude>=0,_ratio>=0,_nb_iter>=0,_blend_operator={_+_|_max_|_min_},_kernel={_0=deriche_|_1=gaussian_|_2=box_|_3=triangle_|_4=quadratic_},_normalize_scales={_0_|_1_},_axes") ); return 0;;
        "blur_linear" | "-blur_linear" | "+blur_linear")
            COMPREPLY=( $(compgen -W "> amplitude1[%],_amplitude2[%],_angle,_boundary_conditions={_0=dirichlet_|_1=neumann_}") ); return 0;;
        "blur_radial" | "-blur_radial" | "+blur_radial")
            COMPREPLY=( $(compgen -W "> amplitude[%],_center_x[%],_center_y[%]") ); return 0;;
        "blur_selective" | "-blur_selective" | "+blur_selective")
            COMPREPLY=( $(compgen -W "> sigma>=0,_edges>0,_nb_scales>0") ); return 0;;
        "blur_x" | "-blur_x" | "+blur_x")
            COMPREPLY=( $(compgen -W "> amplitude[%]>=0,_boundary_conditions={_0=dirichlet_|_1=neumann_|_2=periodic_|_3=mirror_}") ); return 0;;
        "blur_xy" | "-blur_xy" | "+blur_xy")
            COMPREPLY=( $(compgen -W "> amplitude_x[%],amplitude_y[%],_boundary_conditions={_0=dirichlet_|_1=neumann_|_2=periodic_|_3=mirror_}") ); return 0;;
        "blur_xyz" | "-blur_xyz" | "+blur_xyz")
            COMPREPLY=( $(compgen -W "> amplitude_x[%],amplitude_y[%],amplitude_z,_boundary_conditions={_0=dirichlet_|_1=neumann_|_2=periodic_|_3=mirror_}") ); return 0;;
        "blur_y" | "-blur_y" | "+blur_y")
            COMPREPLY=( $(compgen -W "> amplitude[%]>=0,_boundary_conditions={_0=dirichlet_|_1=neumann_|_2=periodic_|_3=mirror_}") ); return 0;;
        "blur_z" | "-blur_z" | "+blur_z")
            COMPREPLY=( $(compgen -W "> amplitude[%]>=0,_boundary_conditions={_0=dirichlet_|_1=neumann_|_2=periodic_|_3=mirror_}") ); return 0;;
        "box3d" | "-box3d" | "+box3d")
            COMPREPLY=( $(compgen -W "> _size_x,_size_y,_size_z") ); return 0;;
        "boxfilter" | "-boxfilter" | "+boxfilter")
            COMPREPLY=( $(compgen -W "size>=0[%],_order,_boundary_conditions,_nb_iter>=0 axes,size>=0[%],_order,_boundary_conditions,_nb_iter>=0") ); return 0;;
        "boxfitting" | "-boxfitting" | "+boxfitting")
            COMPREPLY=( $(compgen -W "> _min_box_size>=1,_max_box_size>=0,_initial_density>=0,_min_spacing>0") ); return 0;;
        "brushify" | "-brushify" | "+brushify")
            COMPREPLY=( $(compgen -W "> [brush],_brush_nb_sizes>=1,0<=_brush_min_size_factor<=1,_brush_nb_orientations>=1,_brush_light_type,0<=_brush_light_strength<=1,_brush_opacity,_painting_density[%]>=0,0<=_painting_contours_coherence<=1,0<=_painting_orientation_coherence<=1,_painting_coherence_alpha[%]>=0,_painting_coherence_sigma[%]>=0,_painting_primary_angle,0<=_painting_angle_dispersion<=1") ); return 0;;
        "bsl" | "-bsl" | "+bsl")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "bsr" | "-bsr" | "+bsr")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "c" | "-c" | "+c")
            COMPREPLY=( $(compgen -W "{_value0[%]_|_[image0]_},{_value1[%]_|_[image1]_} [image]") ); return 0;;
        "camera" | "-camera" | "+camera")
            COMPREPLY=( $(compgen -W "> _camera_index>=0,_nb_frames>0,_skip_frames>=0,_capture_width>=0,_capture_height>=0") ); return 0;;
        "cartoon" | "-cartoon" | "+cartoon")
            COMPREPLY=( $(compgen -W "> _smoothness,_sharpening,_threshold>=0,_thickness>=0,_color>=0,quantization>0") ); return 0;;
        "cast" | "-cast" | "+cast")
            COMPREPLY=( $(compgen -W "> datatype_source,datatype_target") ); return 0;;
        "channels" | "-channels" | "+channels")
            COMPREPLY=( $(compgen -W "> c0[%],_c1[%]") ); return 0;;
        "check3d" | "-check3d" | "+check3d")
            COMPREPLY=( $(compgen -W "> _is_full_check={_0_|_1_}") ); return 0;;
        "check" | "-check" | "+check")
            COMPREPLY=( $(compgen -W "> condition") ); return 0;;
        "chessboard" | "-chessboard" | "+chessboard")
            COMPREPLY=( $(compgen -W "> size1>0,_size2>0,_offset1,_offset2,_angle,_opacity,_color1,...,_color2,...") ); return 0;;
        "circle3d" | "-circle3d" | "+circle3d")
            COMPREPLY=( $(compgen -W "> _x0,_y0,_z0,_radius>=0") ); return 0;;
        "circle" | "-circle" | "+circle")
            COMPREPLY=( $(compgen -W "> x[%],y[%],R[%],_opacity,_pattern,_color1,...") ); return 0;;
        "circles3d" | "-circles3d" | "+circles3d")
            COMPREPLY=( $(compgen -W "> _radius>=0,_is_wireframe={_0_|_1_}") ); return 0;;
        "close_binary" | "-close_binary" | "+close_binary")
            COMPREPLY=( $(compgen -W "> 0<=_endpoint_rate<=100,_endpoint_connectivity>=0,_spline_distmax>=0,_segment_distmax>=0,0<=_spline_anglemax<=180,_spline_roundness>=0,_area_min>=0,_allow_self_intersection={_0_|_1_}") ); return 0;;
        "closing" | "-closing" | "+closing")
            COMPREPLY=( $(compgen -W "size>=0 size_x>=0,size_y>=0,_size_z>=0 [kernel],_boundary_conditions,_is_real={_0=binary-mode_|_1=real-mode_}") ); return 0;;
        "closing_circ" | "-closing_circ" | "+closing_circ")
            COMPREPLY=( $(compgen -W "> _size>=0,_is_real={_0_|_1_}") ); return 0;;
        "clut" | "-clut" | "+clut")
            COMPREPLY=( $(compgen -W "> \"clut_name\",_resolution>0,_cut_and_round={_0=no_|_1=yes_}") ); return 0;;
        "col3d" | "-col3d" | "+col3d")
            COMPREPLY=( $(compgen -W "R,_G,_B,_opacity (no_arg)") ); return 0;;
        "color3d" | "-color3d" | "+color3d")
            COMPREPLY=( $(compgen -W "R,_G,_B,_opacity (no_arg)") ); return 0;;
        "color_ellipses" | "-color_ellipses" | "+color_ellipses")
            COMPREPLY=( $(compgen -W "> _count>0,_radius>=0,_opacity>=0") ); return 0;;
        "colorblind" | "-colorblind" | "+colorblind")
            COMPREPLY=( $(compgen -W "> type={_0=protanopia_|_1=protanomaly_|_2=deuteranopia_|_3=deuteranomaly_|_4=tritanopia_|_5=tritanomaly_|_6=achromatopsia_|_7=achromatomaly_}") ); return 0;;
        "colormap" | "-colormap" | "+colormap")
            COMPREPLY=( $(compgen -W "> nb_levels>=0,_method={_0=median-cut_|_1=k-means_},_sort_vectors") ); return 0;;
        "columns" | "-columns" | "+columns")
            COMPREPLY=( $(compgen -W "> x0[%],_x1[%]") ); return 0;;
        "command" | "-command" | "+command")
            COMPREPLY=( $(compgen -W "_add_debug_info={_0_|_1_},{_filename_|_http[s] //URL_|_\"string\"_}") ); return 0;;
        "compress_clut" | "-compress_clut" | "+compress_clut")
            COMPREPLY=( $(compgen -W "> _max_error>0,_avg_error>0,_max_nbpoints>=8_|_0_(unlimited),_error_metric={_0=L2-norm_|_1=deltaE_1976_|_2=deltaE_2000_},_reconstruction_colorspace={_0=srgb_|_1=rgb_|_2=lab_},_try_rbf_first={_0_|_1_}") ); return 0;;
        "compress_rle" | "-compress_rle" | "+compress_rle")
            COMPREPLY=( $(compgen -W "> _is_binary_data={_0_|_1_},_maximum_sequence_length>=0") ); return 0;;
        "cone3d" | "-cone3d" | "+cone3d")
            COMPREPLY=( $(compgen -W "> _radius,_height,_nb_subdivisions>0") ); return 0;;
        "convolve" | "-convolve" | "+convolve")
            COMPREPLY=( $(compgen -W "> [mask],_boundary_conditions,_is_normalized={_0_|_1_},_channel_mode,_xcenter,_ycenter,_zcenter,_xstart,_ystart,_zstart,_xend,_yend,_zend,_xstride,_ystride,_zstride,_xdilation,_ydilation,_zdilation,interpolation_type") ); return 0;;
        "convolve_fft" | "-convolve_fft" | "+convolve_fft")
            COMPREPLY=( $(compgen -W "> [mask],_boundary_conditions") ); return 0;;
        "correlate" | "-correlate" | "+correlate")
            COMPREPLY=( $(compgen -W "> [mask],_boundary_conditions,_is_normalized={_0_|_1_},_channel_mode,_xcenter,_ycenter,_zcenter,_xstart,_ystart,_zstart,_xend,_yend,_zend,_xstride,_ystride,_zstride,_xdilation,_ydilation,_zdilation,interpolation_type") ); return 0;;
        "covariance_vectors" | "-covariance_vectors" | "+covariance_vectors")
            COMPREPLY=( $(compgen -W "> _avg_outvarname") ); return 0;;
        "cracks" | "-cracks" | "+cracks")
            COMPREPLY=( $(compgen -W "> 0<=_density<=100,_is_relief={_0_|_1_},_opacity,_color1,...") ); return 0;;
        "crop" | "-crop" | "+crop")
            COMPREPLY=( $(compgen -W "x0[%],x1[%],_boundary_conditions x0[%],y0[%],x1[%],y1[%],_boundary_conditions x0[%],y0[%],z0[%],x1[%],y1[%],z1[%],_boundary_conditions x0[%],y0[%],z0[%],c0[%],x1[%],y1[%],z1[%],c1[%],_boundary_conditions") ); return 0;;
        "cross_correlation" | "-cross_correlation" | "+cross_correlation")
            COMPREPLY=( $(compgen -W "> [mask]") ); return 0;;
        "cubes3d" | "-cubes3d" | "+cubes3d")
            COMPREPLY=( $(compgen -W "> _size>=0") ); return 0;;
        "cubism" | "-cubism" | "+cubism")
            COMPREPLY=( $(compgen -W "> _density>=0,0<=_thickness<=50,_max_angle,_opacity,_smoothness>=0") ); return 0;;
        "cumulate" | "-cumulate" | "+cumulate")
            COMPREPLY=( $(compgen -W "{_x_|_y_|_z_|_c_}...{_x_|_y_|_z_|_c_} (no_arg)") ); return 0;;
        "cup3d" | "-cup3d" | "+cup3d")
            COMPREPLY=( $(compgen -W "> _resolution>0") ); return 0;;
        "cursor" | "-cursor" | "+cursor")
            COMPREPLY=( $(compgen -W "> _mode_=_{_0=hide_|_1=show_}") ); return 0;;
        "cut" | "-cut" | "+cut")
            COMPREPLY=( $(compgen -W "{_value0[%]_|_[image0]_},{_value1[%]_|_[image1]_} [image]") ); return 0;;
        "cylinder3d" | "-cylinder3d" | "+cylinder3d")
            COMPREPLY=( $(compgen -W "> _radius,_height,_nb_subdivisions>0") ); return 0;;
        "d3d" | "-d3d" | "+d3d")
            COMPREPLY=( $(compgen -W "_[background_image],_exit_on_anykey={_0_|_1_} _exit_on_anykey={_0_|_1_}") ); return 0;;
        "d" | "-d" | "+d")
            COMPREPLY=( $(compgen -W "> _X[%]>=0,_Y[%]>=0,_Z[%]>=0,_exit_on_anykey={_0_|_1_}") ); return 0;;
        "da" | "-da" | "+da")
            COMPREPLY=( $(compgen -W "> _width>0,_height>0") ); return 0;;
        "db3d" | "-db3d" | "+db3d")
            COMPREPLY=( $(compgen -W "> _is_double_sided={_0_|_1_}") ); return 0;;
        "dct" | "-dct" | "+dct")
            COMPREPLY=( $(compgen -W "_{_x_|_y_|_z_}...{_x_|_y_|_z_} (no_arg)") ); return 0;;
        "deblur" | "-deblur" | "+deblur")
            COMPREPLY=( $(compgen -W "> amplitude[%]>=0,_nb_iter>=0,_dt>=0,_regul>=0,_regul_type={_0=Tikhonov_|_1=meancurv._|_2=TV_}") ); return 0;;
        "deblur_goldmeinel" | "-deblur_goldmeinel" | "+deblur_goldmeinel")
            COMPREPLY=( $(compgen -W "> sigma>=0,_nb_iter>=0,_acceleration>=0,_kernel_type={_0=deriche_|_1=gaussian_}.") ); return 0;;
        "deblur_richardsonlucy" | "-deblur_richardsonlucy" | "+deblur_richardsonlucy")
            COMPREPLY=( $(compgen -W "> sigma>=0,_nb_iter>=0,__kernel_type={_0=deriche_|_1=gaussian_}.") ); return 0;;
        "dec2bin" | "-dec2bin" | "+dec2bin")
            COMPREPLY=( $(compgen -W "> decimal_int1,...") ); return 0;;
        "dec2hex" | "-dec2hex" | "+dec2hex")
            COMPREPLY=( $(compgen -W "> decimal_int1,...") ); return 0;;
        "dec2oct" | "-dec2oct" | "+dec2oct")
            COMPREPLY=( $(compgen -W "> decimal_int1,...") ); return 0;;
        "dec2str" | "-dec2str" | "+dec2str")
            COMPREPLY=( $(compgen -W "> decimal_int1,...") ); return 0;;
        "dec" | "-dec" | "+dec")
            COMPREPLY=( $(compgen -W "> decimal_int1,...") ); return 0;;
        "decompress_clut" | "-decompress_clut" | "+decompress_clut")
            COMPREPLY=( $(compgen -W "> _width>0,_height>0,_depth>0,_reconstruction_colorspace={_0=srgb_|_1=rgb_|_2=lab_}") ); return 0;;
        "decompress_clut_pde" | "-decompress_clut_pde" | "+decompress_clut_pde")
            COMPREPLY=( $(compgen -W "> _width>0,_height>0,_depth>0,_reconstruction_colorspace={_0=srgb_|_1=rgb_|_2=lab_}") ); return 0;;
        "decompress_clut_rbf" | "-decompress_clut_rbf" | "+decompress_clut_rbf")
            COMPREPLY=( $(compgen -W "> _width>0,_height>0,_depth>0,_reconstruction_colorspace={_0=srgb_|_1=rgb_|_2=lab_}") ); return 0;;
        "deconvolve_fft" | "-deconvolve_fft" | "+deconvolve_fft")
            COMPREPLY=( $(compgen -W "> [kernel],_regularization>=0") ); return 0;;
        "deform" | "-deform" | "+deform")
            COMPREPLY=( $(compgen -W "> _amplitude>=0,_interpolation") ); return 0;;
        "deinterlace" | "-deinterlace" | "+deinterlace")
            COMPREPLY=( $(compgen -W "> _method={_0_|_1_}") ); return 0;;
        "delaunay" | "-delaunay" | "+delaunay")
            COMPREPLY=( $(compgen -W "> _output_type={_0=image_|_1=coordinates/triangles_}") ); return 0;;
        "delete" | "-delete" | "+delete")
            COMPREPLY=( $(compgen -W "> filename1[,filename2,...]") ); return 0;;
        "deltaE" | "-deltaE" | "+deltaE")
            COMPREPLY=( $(compgen -W "> [image],_metric={_0=deltaE_1976_|_1=deltaE_2000_},\"_to_Lab_command\"") ); return 0;;
        "demos" | "-demos" | "+demos")
            COMPREPLY=( $(compgen -W "> _run_in_parallel={_0=no_|_1=yes_|_2=auto_}") ); return 0;;
        "denoise" | "-denoise" | "+denoise")
            COMPREPLY=( $(compgen -W "[guide],std_deviation_s[%]>=0,_std_deviation_r[%]>=0,_patch_size>0,_lookup_size>0,_smoothness,_fast_approx={_0_|_1_} std_deviation_s[%]>=0,_std_deviation_r[%]>=0,_patch_size>0,_lookup_size>0,_smoothness,_fast_approx={_0_|_1_}") ); return 0;;
        "denoise_cnn" | "-denoise_cnn" | "+denoise_cnn")
            COMPREPLY=( $(compgen -W "> _noise_type={_0=soft_|_1=heavy_|_2=heavy_(faster)_|_3=poisson+gaussian_|_4=poisson+gaussian2_},_patch_size>0") ); return 0;;
        "denoise_haar" | "-denoise_haar" | "+denoise_haar")
            COMPREPLY=( $(compgen -W "> _threshold>=0,_nb_scales>=0,_cycle_spinning>0") ); return 0;;
        "denoise_patchpca" | "-denoise_patchpca" | "+denoise_patchpca")
            COMPREPLY=( $(compgen -W "> _strength>=0,_patch_size>0,_lookup_size>0,_spatial_sampling>0") ); return 0;;
        "deriche" | "-deriche" | "+deriche")
            COMPREPLY=( $(compgen -W "> std_deviation>=0[%],order={_0_|_1_|_2_},axis={_x_|_y_|_z_|_c_},_boundary_conditions") ); return 0;;
        "detect_skin" | "-detect_skin" | "+detect_skin")
            COMPREPLY=( $(compgen -W "> 0<=tolerance<=1,_skin_x,_skin_y,_skin_radius>=0") ); return 0;;
        "dg" | "-dg" | "+dg")
            COMPREPLY=( $(compgen -W "> _width>=0,_height>=0,_plot_type,_vertex_type,_xmin,_xmax,_ymin,_ymax,_xlabel,_ylabel") ); return 0;;
        "dh" | "-dh" | "+dh")
            COMPREPLY=( $(compgen -W "> _width>=0,_height>=0,_clusters>0,_min_value[%],_max_value[%],_show_axes={_0_|_1_},_expression.") ); return 0;;
        "diffusiontensors" | "-diffusiontensors" | "+diffusiontensors")
            COMPREPLY=( $(compgen -W "> _sharpness>=0,0<=_anisotropy<=1,_alpha[%],_sigma[%],is_sqrt={_0_|_1_}") ); return 0;;
        "dijkstra" | "-dijkstra" | "+dijkstra")
            COMPREPLY=( $(compgen -W "> starting_node>=0,ending_node>=0") ); return 0;;
        "dilate" | "-dilate" | "+dilate")
            COMPREPLY=( $(compgen -W "size>=0 size_x>=0,size_y>=0,size_z>=0 [kernel],_boundary_conditions,_is_real={_0=binary-mode_|_1=real-mode_}") ); return 0;;
        "dilate_circ" | "-dilate_circ" | "+dilate_circ")
            COMPREPLY=( $(compgen -W "> _size>=0,_boundary_conditions,_is_real={_0_|_1_}") ); return 0;;
        "dilate_oct" | "-dilate_oct" | "+dilate_oct")
            COMPREPLY=( $(compgen -W "> _size>=0,_boundary_conditions,_is_real={_0_|_1_}") ); return 0;;
        "dilate_threshold" | "-dilate_threshold" | "+dilate_threshold")
            COMPREPLY=( $(compgen -W "> size_x>=1,size_y>=1,size_z>=1,_threshold>=0,_boundary_conditions") ); return 0;;
        "discard" | "-discard" | "+discard")
            COMPREPLY=( $(compgen -W "_value1,_value2,... {_x_|_y_|_z_|_c}...{_x_|_y_|_z_|_c},_value1,_value2,... (no_arg)") ); return 0;;
        "displacement" | "-displacement" | "+displacement")
            COMPREPLY=( $(compgen -W "> [source_image],_smoothness,_precision>=0,_nb_scales>=0,_iteration_max>=0,is_backward={_0_|_1_},_[guide]") ); return 0;;
        "display3d" | "-display3d" | "+display3d")
            COMPREPLY=( $(compgen -W "_[background_image],_exit_on_anykey={_0_|_1_} _exit_on_anykey={_0_|_1_}") ); return 0;;
        "display" | "-display" | "+display")
            COMPREPLY=( $(compgen -W "> _X[%]>=0,_Y[%]>=0,_Z[%]>=0,_exit_on_anykey={_0_|_1_}") ); return 0;;
        "display_array" | "-display_array" | "+display_array")
            COMPREPLY=( $(compgen -W "> _width>0,_height>0") ); return 0;;
        "display_graph" | "-display_graph" | "+display_graph")
            COMPREPLY=( $(compgen -W "> _width>=0,_height>=0,_plot_type,_vertex_type,_xmin,_xmax,_ymin,_ymax,_xlabel,_ylabel") ); return 0;;
        "display_histogram" | "-display_histogram" | "+display_histogram")
            COMPREPLY=( $(compgen -W "> _width>=0,_height>=0,_clusters>0,_min_value[%],_max_value[%],_show_axes={_0_|_1_},_expression.") ); return 0;;
        "display_parametric" | "-display_parametric" | "+display_parametric")
            COMPREPLY=( $(compgen -W "> _width>0,_height>0,_outline_opacity,_vertex_radius>=0,_is_antialiased={_0_|_1_},_is_decorated={_0_|_1_},_xlabel,_ylabel") ); return 0;;
        "display_polar" | "-display_polar" | "+display_polar")
            COMPREPLY=( $(compgen -W "> _width>32,_height>32,_outline_type,_fill_R,_fill_G,_fill_B,_theta_start,_theta_end,_xlabel,_ylabel") ); return 0;;
        "display_quiver" | "-display_quiver" | "+display_quiver")
            COMPREPLY=( $(compgen -W "> _size_factor>0,_arrow_size>=0,_color_mode={_0=monochrome_|_1=grayscale_|_2=color_}") ); return 0;;
        "display_rgba" | "-display_rgba" | "+display_rgba")
            COMPREPLY=( $(compgen -W "> _background_RGB_color") ); return 0;;
        "display_tensors" | "-display_tensors" | "+display_tensors")
            COMPREPLY=( $(compgen -W "> _size_factor>0,_ellipse_size>=0,_color_mode={_0=monochrome_|_1=grayscale_|_2=color_},_outline>=0") ); return 0;;
        "display_warp" | "-display_warp" | "+display_warp")
            COMPREPLY=( $(compgen -W "> _cell_size>0") ); return 0;;
        "distance" | "-distance" | "+distance")
            COMPREPLY=( $(compgen -W "isovalue[%],_metric isovalue[%],[metric],_method") ); return 0;;
        "div3d" | "-div3d" | "+div3d")
            COMPREPLY=( $(compgen -W "factor factor_x,factor_y,_factor_z") ); return 0;;
        "div" | "-div" | "+div")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "div_complex" | "-div_complex" | "+div_complex")
            COMPREPLY=( $(compgen -W "> [divider_real,divider_imag],_epsilon>=0") ); return 0;;
        "dog" | "-dog" | "+dog")
            COMPREPLY=( $(compgen -W "> _sigma1>=0[%],_sigma2>=0[%]") ); return 0;;
        "double3d" | "-double3d" | "+double3d")
            COMPREPLY=( $(compgen -W "> _is_double_sided={_0_|_1_}") ); return 0;;
        "dq" | "-dq" | "+dq")
            COMPREPLY=( $(compgen -W "> _size_factor>0,_arrow_size>=0,_color_mode={_0=monochrome_|_1=grayscale_|_2=color_}") ); return 0;;
        "draw_whirl" | "-draw_whirl" | "+draw_whirl")
            COMPREPLY=( $(compgen -W "> _amplitude>=0") ); return 0;;
        "drawing" | "-drawing" | "+drawing")
            COMPREPLY=( $(compgen -W "> _amplitude>=0") ); return 0;;
        "drgba" | "-drgba" | "+drgba")
            COMPREPLY=( $(compgen -W "> _background_RGB_color") ); return 0;;
        "drop_shadow" | "-drop_shadow" | "+drop_shadow")
            COMPREPLY=( $(compgen -W "> _offset_x[%],_offset_y[%],_smoothness[%]>=0,0<=_curvature<=1,_expand_size={_0_|_1_}") ); return 0;;
        "dt" | "-dt" | "+dt")
            COMPREPLY=( $(compgen -W "> _size_factor>0,_ellipse_size>=0,_color_mode={_0=monochrome_|_1=grayscale_|_2=color_},_outline>=0") ); return 0;;
        "dw" | "-dw" | "+dw")
            COMPREPLY=( $(compgen -W "> _cell_size>0") ); return 0;;
        "e" | "-e" | "+e")
            COMPREPLY=( $(compgen -W "> message") ); return 0;;
        "echo" | "-echo" | "+echo")
            COMPREPLY=( $(compgen -W "> message") ); return 0;;
        "echo_file" | "-echo_file" | "+echo_file")
            COMPREPLY=( $(compgen -W "> filename,message") ); return 0;;
        "edgels" | "-edgels" | "+edgels")
            COMPREPLY=( $(compgen -W "x0,y0 (no_arg)") ); return 0;;
        "edges" | "-edges" | "+edges")
            COMPREPLY=( $(compgen -W "> _threshold[%]>=0") ); return 0;;
        "elevate" | "-elevate" | "+elevate")
            COMPREPLY=( $(compgen -W "> _depth,_is_plain={_0_|_1_},_is_colored={_0_|_1_}") ); return 0;;
        "elevation3d" | "-elevation3d" | "+elevation3d")
            COMPREPLY=( $(compgen -W "{_z-factor_|_[elevation_map]_|_\'formula\'_},base_height={_-1_|_>=0_} (no_arg)") ); return 0;;
        "elif" | "-elif" | "+elif")
            COMPREPLY=( $(compgen -W "> condition") ); return 0;;
        "ellipse" | "-ellipse" | "+ellipse")
            COMPREPLY=( $(compgen -W "> x[%],y[%],R[%],r[%],_angle,_opacity,_pattern,_color1,...") ); return 0;;
        "ellipsionism" | "-ellipsionism" | "+ellipsionism")
            COMPREPLY=( $(compgen -W "> _R>0[%],_r>0[%],_smoothness>=0[%],_opacity,_outline>0,_density>0") ); return 0;;
        "endian" | "-endian" | "+endian")
            COMPREPLY=( $(compgen -W "> _datatype") ); return 0;;
        "eq" | "-eq" | "+eq")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "equalize" | "-equalize" | "+equalize")
            COMPREPLY=( $(compgen -W "> _nb_levels>0[%],_value_min[%],_value_max[%]") ); return 0;;
        "erode" | "-erode" | "+erode")
            COMPREPLY=( $(compgen -W "size>=0 size_x>=0,size_y>=0,_size_z>=0 [kernel],_boundary_conditions,_is_real={_0=binary-mode_|_1=real-mode_}") ); return 0;;
        "erode_circ" | "-erode_circ" | "+erode_circ")
            COMPREPLY=( $(compgen -W "> _size>=0,_boundary_conditions,_is_real={_0_|_1_}") ); return 0;;
        "erode_oct" | "-erode_oct" | "+erode_oct")
            COMPREPLY=( $(compgen -W "> _size>=0,_boundary_conditions,_is_real={_0_|_1_}") ); return 0;;
        "erode_threshold" | "-erode_threshold" | "+erode_threshold")
            COMPREPLY=( $(compgen -W "> size_x>=1,size_y>=1,size_z>=1,_threshold>=0,_boundary_conditions") ); return 0;;
        "error" | "-error" | "+error")
            COMPREPLY=( $(compgen -W "> message") ); return 0;;
        "euclidean2polar" | "-euclidean2polar" | "+euclidean2polar")
            COMPREPLY=( $(compgen -W "> _center_x[%],_center_y[%],_stretch_factor>0,_boundary_conditions={_0=dirichlet_|_1=neumann_|_2=periodic_|_3=mirror_}") ); return 0;;
        "eval" | "-eval" | "+eval")
            COMPREPLY=( $(compgen -W "> expression") ); return 0;;
        "exec" | "-exec" | "+exec")
            COMPREPLY=( $(compgen -W "> _is_verbose={_0_|_1_},\"command\"") ); return 0;;
        "exec_out" | "-exec_out" | "+exec_out")
            COMPREPLY=( $(compgen -W "> _mode,\"command\"") ); return 0;;
        "expand_x" | "-expand_x" | "+expand_x")
            COMPREPLY=( $(compgen -W "> size_x>=0,_boundary_conditions={_0=dirichlet_|_1=neumann_|_2=periodic_|_3=mirror_}") ); return 0;;
        "expand_xy" | "-expand_xy" | "+expand_xy")
            COMPREPLY=( $(compgen -W "> size>=0,_boundary_conditions={_0=dirichlet_|_1=neumann_|_2=periodic_|_3=mirror_}") ); return 0;;
        "expand_xyz" | "-expand_xyz" | "+expand_xyz")
            COMPREPLY=( $(compgen -W "> size>=0,_boundary_conditions={_0=dirichlet_|_1=neumann_|_2=periodic_|_3=mirror_}") ); return 0;;
        "expand_y" | "-expand_y" | "+expand_y")
            COMPREPLY=( $(compgen -W "> size_y>=0,_boundary_conditions={_0=dirichlet_|_1=neumann_|_2=periodic_|_3=mirror_}") ); return 0;;
        "expand_z" | "-expand_z" | "+expand_z")
            COMPREPLY=( $(compgen -W "> size_z>=0,_boundary_conditions={_0=dirichlet_|_1=neumann_|_2=periodic_|_3=mirror_}") ); return 0;;
        "extract" | "-extract" | "+extract")
            COMPREPLY=( $(compgen -W "> \"condition\",_output_type={_0=xyzc-coordinates_|_1=xyz-coordinates_|_2=scalar-values_|_3=vector-values_}") ); return 0;;
        "extract_region" | "-extract_region" | "+extract_region")
            COMPREPLY=( $(compgen -W "> [label_image],_extract_xyz_coordinates={_0_|_1_},_label_1,...,_label_M") ); return 0;;
        "extrude3d" | "-extrude3d" | "+extrude3d")
            COMPREPLY=( $(compgen -W "> _depth>0,_resolution>0,_smoothness[%]>=0") ); return 0;;
        "eye" | "-eye" | "+eye")
            COMPREPLY=( $(compgen -W "> _size>0") ); return 0;;
        "f3d" | "-f3d" | "+f3d")
            COMPREPLY=( $(compgen -W "> focale") ); return 0;;
        "f" | "-f" | "+f")
            COMPREPLY=( $(compgen -W "value1,_value2,... [image] \'formula\'") ); return 0;;
        "fact" | "-fact" | "+fact")
            COMPREPLY=( $(compgen -W "> value") ); return 0;;
        "fade_diamond" | "-fade_diamond" | "+fade_diamond")
            COMPREPLY=( $(compgen -W "> 0<=_start<=100,0<=_end<=100") ); return 0;;
        "fade_files" | "-fade_files" | "+fade_files")
            COMPREPLY=( $(compgen -W "> \"filename_pattern\",_nb_inner_frames>0,_first_frame>=0,_last_frame={_>=0_|_-1=last_},_frame_step>=1,_output_filename") ); return 0;;
        "fade_linear" | "-fade_linear" | "+fade_linear")
            COMPREPLY=( $(compgen -W "> _angle,0<=_start<=100,0<=_end<=100") ); return 0;;
        "fade_radial" | "-fade_radial" | "+fade_radial")
            COMPREPLY=( $(compgen -W "> 0<=_start<=100,0<=_end<=100") ); return 0;;
        "fade_video" | "-fade_video" | "+fade_video")
            COMPREPLY=( $(compgen -W "> video_filename,_nb_inner_frames>0,_first_frame>=0,_last_frame={_>=0_|_-1=last_},_frame_step>=1,_output_filename") ); return 0;;
        "fade_x" | "-fade_x" | "+fade_x")
            COMPREPLY=( $(compgen -W "> 0<=_start<=100,0<=_end<=100") ); return 0;;
        "fade_y" | "-fade_y" | "+fade_y")
            COMPREPLY=( $(compgen -W "> 0<=_start<=100,0<=_end<=100") ); return 0;;
        "fade_z" | "-fade_z" | "+fade_z")
            COMPREPLY=( $(compgen -W "> 0<=_start<=100,0<=_end<=100") ); return 0;;
        "fc" | "-fc" | "+fc")
            COMPREPLY=( $(compgen -W "> col1,...,colN") ); return 0;;
        "fft" | "-fft" | "+fft")
            COMPREPLY=( $(compgen -W "> _{_x_|_y_|_z_}...{_x_|_y_|_z_}") ); return 0;;
        "fibonacci" | "-fibonacci" | "+fibonacci")
            COMPREPLY=( $(compgen -W "> N>=0") ); return 0;;
        "file_mv" | "-file_mv" | "+file_mv")
            COMPREPLY=( $(compgen -W "> filename_src,filename_dest") ); return 0;;
        "filename" | "-filename" | "+filename")
            COMPREPLY=( $(compgen -W "> filename,_number1,_number2,...,_numberN") ); return 0;;
        "files2img" | "-files2img" | "+files2img")
            COMPREPLY=( $(compgen -W "> _mode,path") ); return 0;;
        "files2video" | "-files2video" | "+files2video")
            COMPREPLY=( $(compgen -W "> \"filename_pattern\",_output_filename,_fps>0,_codec") ); return 0;;
        "files" | "-files" | "+files")
            COMPREPLY=( $(compgen -W "> _mode,path") ); return 0;;
        "fill" | "-fill" | "+fill")
            COMPREPLY=( $(compgen -W "value1,_value2,... [image] \'formula\'") ); return 0;;
        "fill_color" | "-fill_color" | "+fill_color")
            COMPREPLY=( $(compgen -W "> col1,...,colN") ); return 0;;
        "fire_edges" | "-fire_edges" | "+fire_edges")
            COMPREPLY=( $(compgen -W "> _edges>=0,0<=_attenuation<=1,_smoothness>=0,_threshold>=0,_nb_frames>0,_starting_frame>=0,frame_skip>=0") ); return 0;;
        "fisheye" | "-fisheye" | "+fisheye")
            COMPREPLY=( $(compgen -W "> _center_x,_center_y,0<=_radius<=100,_amplitude>=0") ); return 0;;
        "fitratio_wh" | "-fitratio_wh" | "+fitratio_wh")
            COMPREPLY=( $(compgen -W "> min_width,min_height,ratio_wh") ); return 0;;
        "fitscreen" | "-fitscreen" | "+fitscreen")
            COMPREPLY=( $(compgen -W "width,height,_depth,_minimal_size[%],_maximal_size[%] [image],_minimal_size[%],_maximal_size[%]") ); return 0;;
        "flood" | "-flood" | "+flood")
            COMPREPLY=( $(compgen -W "> x[%],_y[%],_z[%],_tolerance>=0,_is_high_connectivity={_0_|_1_},_opacity,_color1,...") ); return 0;;
        "flower" | "-flower" | "+flower")
            COMPREPLY=( $(compgen -W "> _amplitude,_frequency,_offset_r[%],_angle,_center_x[%],_center_y[%],_boundary_conditions={_0=dirichlet_|_1=neumann_|_2=periodic_|_3=mirror}") ); return 0;;
        "focale3d" | "-focale3d" | "+focale3d")
            COMPREPLY=( $(compgen -W "> focale") ); return 0;;
        "for" | "-for" | "+for")
            COMPREPLY=( $(compgen -W "> condition") ); return 0;;
        "fractalize" | "-fractalize" | "+fractalize")
            COMPREPLY=( $(compgen -W "> 0<=detail_level<=1") ); return 0;;
        "frame" | "-frame" | "+frame")
            COMPREPLY=( $(compgen -W "> size_x[%],_size_y[%],_col1,...,_colN") ); return 0;;
        "frame_blur" | "-frame_blur" | "+frame_blur")
            COMPREPLY=( $(compgen -W "> _sharpness>0,_size>=0,_smoothness,_shading,_blur") ); return 0;;
        "frame_cube" | "-frame_cube" | "+frame_cube")
            COMPREPLY=( $(compgen -W "> _depth>=0,_centering_x,_centering_y,_left_side={0=normal_|_1=mirror-x_|_2=mirror-y_|_3=mirror-xy},_right_side,_lower_side,_upper_side") ); return 0;;
        "frame_fuzzy" | "-frame_fuzzy" | "+frame_fuzzy")
            COMPREPLY=( $(compgen -W "> size_x[%]>=0,_size_y[%]>=0,_fuzzyness>=0,_smoothness[%]>=0,_R,_G,_B,_A") ); return 0;;
        "frame_painting" | "-frame_painting" | "+frame_painting")
            COMPREPLY=( $(compgen -W "> _size[%]>=0,0<=_contrast<=1,_profile_smoothness[%]>=0,_R,_G,_B,_vignette_size[%]>=0,_vignette_contrast>=0,_defects_contrast>=0,0<=_defects_density<=100,_defects_size>=0,_defects_smoothness[%]>=0,_serial_number") ); return 0;;
        "frame_pattern" | "-frame_pattern" | "+frame_pattern")
            COMPREPLY=( $(compgen -W "M>=3,_constrain_size={_0_|_1_} M>=3,_[frame_image],_constrain_size={_0_|_1_}") ); return 0;;
        "frame_round" | "-frame_round" | "+frame_round")
            COMPREPLY=( $(compgen -W "> frame_size[%]>=0,radius[%]>=0,_smoothness[%]>=0,_col1,...,_colN") ); return 0;;
        "frame_seamless" | "-frame_seamless" | "+frame_seamless")
            COMPREPLY=( $(compgen -W "> frame_size>=0,_patch_size>0,_blend_size>=0,_frame_direction={_0=inner_(preserve_image_size)_|_1=outer_}") ); return 0;;
        "frame_x" | "-frame_x" | "+frame_x")
            COMPREPLY=( $(compgen -W "> size_x[%],_col1,...,_colN") ); return 0;;
        "frame_xy" | "-frame_xy" | "+frame_xy")
            COMPREPLY=( $(compgen -W "> size_x[%],_size_y[%],_col1,...,_colN") ); return 0;;
        "frame_xyz" | "-frame_xyz" | "+frame_xyz")
            COMPREPLY=( $(compgen -W "> size_x[%],_size_y[%],_size_z[%]_col1,...,_colN") ); return 0;;
        "frame_y" | "-frame_y" | "+frame_y")
            COMPREPLY=( $(compgen -W "> size_y[%],_col1,...,_colN") ); return 0;;
        "function1d" | "-function1d" | "+function1d")
            COMPREPLY=( $(compgen -W "> 0<=smoothness<=1,x0>=0,y0,x1>=0,y1,...,xn>=0,yn") ); return 0;;
        "g" | "-g" | "+g")
            COMPREPLY=( $(compgen -W "{_x_|_y_|_z_|_c_}...{_x_|_y_|_z_|_c_},_scheme,_boundary_conditions (no_arg)") ); return 0;;
        "gaussian" | "-gaussian" | "+gaussian")
            COMPREPLY=( $(compgen -W "> _sigma1[%],_sigma2[%],_angle") ); return 0;;
        "gaussians3d" | "-gaussians3d" | "+gaussians3d")
            COMPREPLY=( $(compgen -W "> _size>0,_opacity") ); return 0;;
        "gcd" | "-gcd" | "+gcd")
            COMPREPLY=( $(compgen -W "> a,b") ); return 0;;
        "ge" | "-ge" | "+ge")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "glow" | "-glow" | "+glow")
            COMPREPLY=( $(compgen -W "> _amplitude>=0") ); return 0;;
        "gmd2ascii" | "-gmd2ascii" | "+gmd2ascii")
            COMPREPLY=( $(compgen -W "_max_line_length>0,_indent_forced_newlines>=0 (no_arg)") ); return 0;;
        "gmd2html" | "-gmd2html" | "+gmd2html")
            COMPREPLY=( $(compgen -W "_include_default_header_footer={_0=none_|_1=Reference_|_2=Tutorial_|_3=News_} (no_arg)") ); return 0;;
        "gradient2rgb" | "-gradient2rgb" | "+gradient2rgb")
            COMPREPLY=( $(compgen -W "> _is_orientation={_0_|_1_}") ); return 0;;
        "gradient" | "-gradient" | "+gradient")
            COMPREPLY=( $(compgen -W "{_x_|_y_|_z_|_c_}...{_x_|_y_|_z_|_c_},_scheme,_boundary_conditions (no_arg)") ); return 0;;
        "gradient_orientation" | "-gradient_orientation" | "+gradient_orientation")
            COMPREPLY=( $(compgen -W "> _dimension={_1_|_2_|_3_}") ); return 0;;
        "graph" | "-graph" | "+graph")
            COMPREPLY=( $(compgen -W "[function_image],_plot_type,_vertex_type,_ymin,_ymax,_opacity,_pattern,_color1,... \'formula\',_resolution>=0,_plot_type,_vertex_type,_xmin,xmax,_ymin,_ymax,_opacity,_pattern,_color1,...") ); return 0;;
        "grid" | "-grid" | "+grid")
            COMPREPLY=( $(compgen -W "> size_x[%]>=0,size_y[%]>=0,_offset_x[%],_offset_y[%],_opacity,_pattern,_color1,...") ); return 0;;
        "gt" | "-gt" | "+gt")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "guided" | "-guided" | "+guided")
            COMPREPLY=( $(compgen -W "[guide],radius[%]>=0,regularization[%]>=0 radius[%]>=0,regularization[%]>=0") ); return 0;;
        "gyroid3d" | "-gyroid3d" | "+gyroid3d")
            COMPREPLY=( $(compgen -W "> _resolution>0,_zoom") ); return 0;;
        "h" | "-h" | "+h")
            COMPREPLY=( $(compgen -W "$coms" -- "$cur") ); return 0;;
        "haar" | "-haar" | "+haar")
            COMPREPLY=( $(compgen -W "> scale>0") ); return 0;;
        "halftone" | "-halftone" | "+halftone")
            COMPREPLY=( $(compgen -W "> nb_levels>=2,_size_dark>=2,_size_bright>=2,_shape={_0=square_|_1=diamond_|_2=circle_|_3=inv-square_|_4=inv-diamond_|_5=inv-circle_},_smoothness[%]>=0") ); return 0;;
        "hardsketchbw" | "-hardsketchbw" | "+hardsketchbw")
            COMPREPLY=( $(compgen -W "> _amplitude>=0,_density>=0,_opacity,0<=_edge_threshold<=100,_is_fast={_0_|_1_}") ); return 0;;
        "hearts" | "-hearts" | "+hearts")
            COMPREPLY=( $(compgen -W "> _density>=0") ); return 0;;
        "heat_flow" | "-heat_flow" | "+heat_flow")
            COMPREPLY=( $(compgen -W "> _nb_iter>=0,_dt,_keep_sequence={_0_|_1_}") ); return 0;;
        "help" | "-help" | "+help")
            COMPREPLY=( $(compgen -W "$coms" -- "$cur") ); return 0;;
        "hessian" | "-hessian" | "+hessian")
            COMPREPLY=( $(compgen -W "{_xx_|_xy_|_xz_|_yy_|_yz_|_zz_}...{_xx_|_xy_|_xz_|_yy_|_yz_|_zz_},_boundary_conditions (no_arg)") ); return 0;;
        "hex2dec" | "-hex2dec" | "+hex2dec")
            COMPREPLY=( $(compgen -W "> hexadecimal_int1,...") ); return 0;;
        "hex2img" | "-hex2img" | "+hex2img")
            COMPREPLY=( $(compgen -W "> \"hexadecimal_string\"") ); return 0;;
        "hex2str" | "-hex2str" | "+hex2str")
            COMPREPLY=( $(compgen -W "> hexadecimal_string") ); return 0;;
        "hex" | "-hex" | "+hex")
            COMPREPLY=( $(compgen -W "> hexadecimal_int1,...") ); return 0;;
        "histogram" | "-histogram" | "+histogram")
            COMPREPLY=( $(compgen -W "> nb_levels>0[%],_min_value[%],_max_value[%]") ); return 0;;
        "histogram_cumul" | "-histogram_cumul" | "+histogram_cumul")
            COMPREPLY=( $(compgen -W "> _nb_levels>0,_is_normalized={_0_|_1_},_val0[%],_val1[%]") ); return 0;;
        "histogram_nd" | "-histogram_nd" | "+histogram_nd")
            COMPREPLY=( $(compgen -W "> nb_levels>0[%],_value0[%],_value1[%]") ); return 0;;
        "histogram_pointwise" | "-histogram_pointwise" | "+histogram_pointwise")
            COMPREPLY=( $(compgen -W "> nb_levels>0[%],_value0[%],_value1[%]") ); return 0;;
        "hough" | "-hough" | "+hough")
            COMPREPLY=( $(compgen -W "> _width>0,_height>0,gradient_norm_voting={_0_|_1_}") ); return 0;;
        "houghsketchbw" | "-houghsketchbw" | "+houghsketchbw")
            COMPREPLY=( $(compgen -W "> _density>=0,_radius>0,0<=_threshold<=100,0<=_opacity<=1,_votesize[%]>0") ); return 0;;
        "idct" | "-idct" | "+idct")
            COMPREPLY=( $(compgen -W "_{_x_|_y_|_z_}...{_x_|_y_|_z_} (no_arg)") ); return 0;;
        "identity" | "-identity" | "+identity")
            COMPREPLY=( $(compgen -W "> _width>=0,_height>=0,_depth>=0") ); return 0;;
        "if" | "-if" | "+if")
            COMPREPLY=( $(compgen -W "> condition") ); return 0;;
        "ifft" | "-ifft" | "+ifft")
            COMPREPLY=( $(compgen -W "> _{_x_|_y_|_z_}...{_x_|_y_|_z_}") ); return 0;;
        "ig" | "-ig" | "+ig")
            COMPREPLY=( $(compgen -W "> pattern") ); return 0;;
        "ihaar" | "-ihaar" | "+ihaar")
            COMPREPLY=( $(compgen -W "> scale>0") ); return 0;;
        "ilaplacian" | "-ilaplacian" | "+ilaplacian")
            COMPREPLY=( $(compgen -W "> {_nb_iterations>0_|_0_},_[initial_estimate]") ); return 0;;
        "image" | "-image" | "+image")
            COMPREPLY=( $(compgen -W "> [sprite],_x[%|~],_y[%|~],_z[%|~],_c[%|~],_opacity,_[opacity_mask],_max_opacity_mask") ); return 0;;
        "imageblocks3d" | "-imageblocks3d" | "+imageblocks3d")
            COMPREPLY=( $(compgen -W "> _maximum_elevation,_smoothness[%]>=0") ); return 0;;
        "imagegrid" | "-imagegrid" | "+imagegrid")
            COMPREPLY=( $(compgen -W "> M>0,_N>0") ); return 0;;
        "imagegrid_hexagonal" | "-imagegrid_hexagonal" | "+imagegrid_hexagonal")
            COMPREPLY=( $(compgen -W "> _resolution>0,0<=_outline<=1") ); return 0;;
        "imagegrid_triangular" | "-imagegrid_triangular" | "+imagegrid_triangular")
            COMPREPLY=( $(compgen -W "> pattern_width>=1,_pattern_height>=1,_pattern_type,0<=_outline_opacity<=1,_outline_color1,...") ); return 0;;
        "imagerubik3d" | "-imagerubik3d" | "+imagerubik3d")
            COMPREPLY=( $(compgen -W "> _xy_tiles>=1,0<=xy_shift<=100,0<=z_shift<=100") ); return 0;;
        "imagesphere3d" | "-imagesphere3d" | "+imagesphere3d")
            COMPREPLY=( $(compgen -W "> _resolution1>=3,_resolution2>=3") ); return 0;;
        "img2ascii" | "-img2ascii" | "+img2ascii")
            COMPREPLY=( $(compgen -W "> _charset,_analysis_scale>0,_analysis_smoothness[%]>=0,_synthesis_scale>0,_output_ascii_filename") ); return 0;;
        "img2base64" | "-img2base64" | "+img2base64")
            COMPREPLY=( $(compgen -W "> _encoding={_0=base64_|_1=base64url_},_store_names={_0_|_1_}") ); return 0;;
        "img2patches" | "-img2patches" | "+img2patches")
            COMPREPLY=( $(compgen -W "> patch_size>0,_overlap[%]>0,_boundary_conditions") ); return 0;;
        "img2text" | "-img2text" | "+img2text")
            COMPREPLY=( $(compgen -W "> _line_separator") ); return 0;;
        "index" | "-index" | "+index")
            COMPREPLY=( $(compgen -W "> {_[palette]_|_palette_name_},0<=_dithering<=1,_map_palette={_0_|_1_}") ); return 0;;
        "inpaint" | "-inpaint" | "+inpaint")
            COMPREPLY=( $(compgen -W "[mask] [mask],0,_fast_method [mask],_patch_size>=1,_lookup_size>=1,_lookup_factor>=0,_lookup_increment!=0,_blend_size>=0,0<=_blend_threshold<=1,_blend_decay>=0,_blend_scales>=1,_is_blend_outer={_0_|_1_}") ); return 0;;
        "inpaint_flow" | "-inpaint_flow" | "+inpaint_flow")
            COMPREPLY=( $(compgen -W "> [mask],_nb_global_iter>=0,_nb_local_iter>=0,_dt>0,_alpha>=0,_sigma>=0") ); return 0;;
        "inpaint_holes" | "-inpaint_holes" | "+inpaint_holes")
            COMPREPLY=( $(compgen -W "> maximal_area[%]>=0,_tolerance>=0,_is_high_connectivity={_0_|_1_}") ); return 0;;
        "inpaint_matchpatch" | "-inpaint_matchpatch" | "+inpaint_matchpatch")
            COMPREPLY=( $(compgen -W "> [mask],_nb_scales={_0=auto_|_>0_},_patch_size>0,_nb_iterations_per_scale>0,_blend_size>=0,_allow_outer_blending={_0_|_1_},_is_already_initialized={_0_|_1_}") ); return 0;;
        "inpaint_morpho" | "-inpaint_morpho" | "+inpaint_morpho")
            COMPREPLY=( $(compgen -W "> [mask]") ); return 0;;
        "inpaint_pde" | "-inpaint_pde" | "+inpaint_pde")
            COMPREPLY=( $(compgen -W "> [mask],_nb_scales[%]>=0,_diffusion_type={_0=isotropic_|_1=Delaunay-guided_|_2=edge-guided_|_3=mask-guided_},_diffusion_iter>=0") ); return 0;;
        "inrange" | "-inrange" | "+inrange")
            COMPREPLY=( $(compgen -W "> min[%],max[%],_include_min_boundary={_0=no_|_1=yes_},_include_max_boundary={_0=no_|_1=yes_}") ); return 0;;
        "invert" | "-invert" | "+invert")
            COMPREPLY=( $(compgen -W "> solver={_0=SVD_|_1=LU_}") ); return 0;;
        "ir" | "-ir" | "+ir")
            COMPREPLY=( $(compgen -W "> min[%],max[%],_include_min_boundary={_0=no_|_1=yes_},_include_max_boundary={_0=no_|_1=yes_}") ); return 0;;
        "is_change" | "-is_change" | "+is_change")
            COMPREPLY=( $(compgen -W "> _value={_0=false_|_1=true_}") ); return 0;;
        "is_ext" | "-is_ext" | "+is_ext")
            COMPREPLY=( $(compgen -W "> filename,_extension") ); return 0;;
        "is_image_arg" | "-is_image_arg" | "+is_image_arg")
            COMPREPLY=( $(compgen -W "> string") ); return 0;;
        "is_pattern" | "-is_pattern" | "+is_pattern")
            COMPREPLY=( $(compgen -W "> string") ); return 0;;
        "is_percent" | "-is_percent" | "+is_percent")
            COMPREPLY=( $(compgen -W "> string") ); return 0;;
        "isoline3d" | "-isoline3d" | "+isoline3d")
            COMPREPLY=( $(compgen -W "isovalue[%] \'formula\',value,_x0,_y0,_x1,_y1,_size_x>0[%],_size_y>0[%]") ); return 0;;
        "isophotes" | "-isophotes" | "+isophotes")
            COMPREPLY=( $(compgen -W "> _nb_levels>0") ); return 0;;
        "isosurface3d" | "-isosurface3d" | "+isosurface3d")
            COMPREPLY=( $(compgen -W "isovalue[%] \'formula\',value,_x0,_y0,_z0,_x1,_y1,_z1,_size_x>0[%],_size_y>0[%],_size_z>0[%]") ); return 0;;
        "j3d" | "-j3d" | "+j3d")
            COMPREPLY=( $(compgen -W "> [object3d],_x[%],_y[%],_z,_opacity,_rendering_mode,_is_double_sided={_0_|_1_},_is_zbuffer={_0_|_1_},_focale,_light_x,_light_y,_light_z,_specular_lightness,_specular_shininess") ); return 0;;
        "j" | "-j" | "+j")
            COMPREPLY=( $(compgen -W "> [sprite],_x[%|~],_y[%|~],_z[%|~],_c[%|~],_opacity,_[opacity_mask],_max_opacity_mask") ); return 0;;
        "jzazbz2rgb" | "-jzazbz2rgb" | "+jzazbz2rgb")
            COMPREPLY=( $(compgen -W "illuminant={_0=D50_|_1=D65_|_2=E_} (no_arg)") ); return 0;;
        "kaleidoscope" | "-kaleidoscope" | "+kaleidoscope")
            COMPREPLY=( $(compgen -W "> _center_x[%],_center_y[%],_radius,_angle,_boundary_conditions={_0=dirichlet_|_1=neumann_|_2=periodic_|_3=mirror_}") ); return 0;;
        "keep_named" | "-keep_named" | "+keep_named")
            COMPREPLY=( $(compgen -W "> \"name1\",\"name2\",...") ); return 0;;
        "kn" | "-kn" | "+kn")
            COMPREPLY=( $(compgen -W "> \"name1\",\"name2\",...") ); return 0;;
        "kuwahara" | "-kuwahara" | "+kuwahara")
            COMPREPLY=( $(compgen -W "> size>0") ); return 0;;
        "l3d" | "-l3d" | "+l3d")
            COMPREPLY=( $(compgen -W "position_x,position_y,position_z [texture] (no_arg)") ); return 0;;
        "lab2rgb" | "-lab2rgb" | "+lab2rgb")
            COMPREPLY=( $(compgen -W "illuminant={_0=D50_|_1=D65_|_2=E_} (no_arg)") ); return 0;;
        "lab2srgb" | "-lab2srgb" | "+lab2srgb")
            COMPREPLY=( $(compgen -W "illuminant={_0=D50_|_1=D65_|_2=E_} (no_arg)") ); return 0;;
        "lab2xyz" | "-lab2xyz" | "+lab2xyz")
            COMPREPLY=( $(compgen -W "illuminant={_0=D50_|_1=D65_|_2=E_} (no_arg)") ); return 0;;
        "lab82rgb" | "-lab82rgb" | "+lab82rgb")
            COMPREPLY=( $(compgen -W "illuminant={_0=D50_|_1=D65_|_2=E_} (no_arg)") ); return 0;;
        "lab82srgb" | "-lab82srgb" | "+lab82srgb")
            COMPREPLY=( $(compgen -W "illuminant={_0=D50_|_1=D65_|_2=E_} (no_arg)") ); return 0;;
        "label3d" | "-label3d" | "+label3d")
            COMPREPLY=( $(compgen -W "> \"text\",font_height>=0,_opacity,_color1,...") ); return 0;;
        "label" | "-label" | "+label")
            COMPREPLY=( $(compgen -W "> _tolerance>=0,is_high_connectivity={_0_|_1_},_is_L2_norm={_0_|_1_}") ); return 0;;
        "label_fg" | "-label_fg" | "+label_fg")
            COMPREPLY=( $(compgen -W "> tolerance>=0,is_high_connectivity={_0_|_1_}") ); return 0;;
        "label_points3d" | "-label_points3d" | "+label_points3d")
            COMPREPLY=( $(compgen -W "> _label_size>0,_opacity") ); return 0;;
        "lathe3d" | "-lathe3d" | "+lathe3d")
            COMPREPLY=( $(compgen -W "> _resolution>0,_smoothness[%]>=0,_max_angle>=0") ); return 0;;
        "lch2rgb" | "-lch2rgb" | "+lch2rgb")
            COMPREPLY=( $(compgen -W "illuminant={_0=D50_|_1=D65_|_2=E_} (no_arg)") ); return 0;;
        "lch82rgb" | "-lch82rgb" | "+lch82rgb")
            COMPREPLY=( $(compgen -W "illuminant={_0=D50_|_1=D65_|_2=E_} (no_arg)") ); return 0;;
        "le" | "-le" | "+le")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "lic" | "-lic" | "+lic")
            COMPREPLY=( $(compgen -W "> _amplitude>0,_channels>0") ); return 0;;
        "light3d" | "-light3d" | "+light3d")
            COMPREPLY=( $(compgen -W "position_x,position_y,position_z [texture] (no_arg)") ); return 0;;
        "light_patch" | "-light_patch" | "+light_patch")
            COMPREPLY=( $(compgen -W "> _density>0,_darkness>=0,_lightness>=0") ); return 0;;
        "light_relief" | "-light_relief" | "+light_relief")
            COMPREPLY=( $(compgen -W "> _ambient_light,_specular_lightness,_specular_size,_darkness,_light_smoothness,_xl,_yl,_zl,_zscale,_opacity_is_heightmap={_0_|_1_}") ); return 0;;
        "lightrays" | "-lightrays" | "+lightrays")
            COMPREPLY=( $(compgen -W "> 100<=_density<=0,_center_x[%],_center_y[%],_ray_length>=0,_ray_attenuation>=0") ); return 0;;
        "line3d" | "-line3d" | "+line3d")
            COMPREPLY=( $(compgen -W "> x0,y0,z0,x1,y1,z1") ); return 0;;
        "line" | "-line" | "+line")
            COMPREPLY=( $(compgen -W "> x0[%],y0[%],x1[%],y1[%],_opacity,_pattern,_color1,...") ); return 0;;
        "linearize_tiles" | "-linearize_tiles" | "+linearize_tiles")
            COMPREPLY=( $(compgen -W "> M>0,_N>0") ); return 0;;
        "linify" | "-linify" | "+linify")
            COMPREPLY=( $(compgen -W "> 0<=_density<=100,_spreading>=0,_resolution[%]>0,_line_opacity>=0,_line_precision>0,_mode={_0=subtractive_|_1=additive_}") ); return 0;;
        "lissajous3d" | "-lissajous3d" | "+lissajous3d")
            COMPREPLY=( $(compgen -W "> resolution>1,a,A,b,B,c,C") ); return 0;;
        "lorem" | "-lorem" | "+lorem")
            COMPREPLY=( $(compgen -W "> _width>0,_height>0") ); return 0;;
        "lt" | "-lt" | "+lt")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "lut_contrast" | "-lut_contrast" | "+lut_contrast")
            COMPREPLY=( $(compgen -W "> _nb_colors>1,_min_rgb_value") ); return 0;;
        "m*" | "-m*" | "+m*")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "m/" | "-m/" | "+m/")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "m3d" | "-m3d" | "+m3d")
            COMPREPLY=( $(compgen -W "> _mode") ); return 0;;
        "mandelbrot" | "-mandelbrot" | "+mandelbrot")
            COMPREPLY=( $(compgen -W "> z0r,z0i,z1r,z1i,_iteration_max>=0,_is_julia={_0_|_1_},_c0r,_c0i,_opacity") ); return 0;;
        "map" | "-map" | "+map")
            COMPREPLY=( $(compgen -W "[palette],_boundary_conditions palette_name,_boundary_conditions") ); return 0;;
        "map_clut" | "-map_clut" | "+map_clut")
            COMPREPLY=( $(compgen -W "> [clut]_|_\"clut_name\"") ); return 0;;
        "map_sphere" | "-map_sphere" | "+map_sphere")
            COMPREPLY=( $(compgen -W "> _width>0,_height>0,_radius,_dilation>0,_fading>=0,_fading_power>=0") ); return 0;;
        "map_sprites" | "-map_sprites" | "+map_sprites")
            COMPREPLY=( $(compgen -W "> _nb_sprites>=1,_allow_rotation={_0=none_|_1=90_deg._|_2=180_deg._}") ); return 0;;
        "map_tones" | "-map_tones" | "+map_tones")
            COMPREPLY=( $(compgen -W "> _threshold>=0,_gamma>=0,_smoothness>=0,nb_iter>=0") ); return 0;;
        "map_tones_fast" | "-map_tones_fast" | "+map_tones_fast")
            COMPREPLY=( $(compgen -W "> _radius[%]>=0,_power>=0") ); return 0;;
        "marble" | "-marble" | "+marble")
            COMPREPLY=( $(compgen -W "> _image_weight,_pattern_weight,_angle,_amplitude,_sharpness>=0,_anisotropy>=0,_alpha,_sigma,_cut_low>=0,_cut_high>=0") ); return 0;;
        "matchpatch" | "-matchpatch" | "+matchpatch")
            COMPREPLY=( $(compgen -W "> [patch_image],patch_width>=1,_patch_height>=1,_patch_depth>=1,_nb_iterations>=0,_nb_randoms>=0,_patch_penalization,_output_score={_0_|_1_},_[guide]") ); return 0;;
        "max" | "-max" | "+max")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "max_patch" | "-max_patch" | "+max_patch")
            COMPREPLY=( $(compgen -W "> _patch_size>=1") ); return 0;;
        "maxabs" | "-maxabs" | "+maxabs")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "maze" | "-maze" | "+maze")
            COMPREPLY=( $(compgen -W "> _width>0,_height>0,_cell_size>0") ); return 0;;
        "maze_mask" | "-maze_mask" | "+maze_mask")
            COMPREPLY=( $(compgen -W "> _cellsize>0") ); return 0;;
        "md3d" | "-md3d" | "+md3d")
            COMPREPLY=( $(compgen -W "> _mode") ); return 0;;
        "mdiv" | "-mdiv" | "+mdiv")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "meancurvature_flow" | "-meancurvature_flow" | "+meancurvature_flow")
            COMPREPLY=( $(compgen -W "> _nb_iter>=0,_dt,_keep_sequence={_0_|_1_}") ); return 0;;
        "median" | "-median" | "+median")
            COMPREPLY=( $(compgen -W "> size>=0,_threshold>0") ); return 0;;
        "median_files" | "-median_files" | "+median_files")
            COMPREPLY=( $(compgen -W "> \"filename_pattern\",_first_frame>=0,_last_frame={_>=0_|_-1=last_},_frame_step>=1,_frame_rows[%]>=1,_is_fast_approximation={_0_|_1_}") ); return 0;;
        "median_video" | "-median_video" | "+median_video")
            COMPREPLY=( $(compgen -W "> video_filename,_first_frame>=0,_last_frame={_>=0_|_-1=last_},_frame_step>=1,_frame_rows[%]>=1,_is_fast_approximation={_0_|_1_}") ); return 0;;
        "meigen" | "-meigen" | "+meigen")
            COMPREPLY=( $(compgen -W "> m>=1") ); return 0;;
        "min" | "-min" | "+min")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "min_patch" | "-min_patch" | "+min_patch")
            COMPREPLY=( $(compgen -W "> _patch_size>=1") ); return 0;;
        "minabs" | "-minabs" | "+minabs")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "minimal_path" | "-minimal_path" | "+minimal_path")
            COMPREPLY=( $(compgen -W "> x0[%]>=0,y0[%]>=0,z0[%]>=0,x1[%]>=0,y1[%]>=0,z1[%]>=0,_is_high_connectivity={_0_|_1_}") ); return 0;;
        "mirror" | "-mirror" | "+mirror")
            COMPREPLY=( $(compgen -W "> {_x_|_y_|_z_}...{_x_|_y_|_z_}") ); return 0;;
        "mix_channels" | "-mix_channels" | "+mix_channels")
            COMPREPLY=( $(compgen -W "(a00,...,aMN) [matrix]") ); return 0;;
        "mix_rgb" | "-mix_rgb" | "+mix_rgb")
            COMPREPLY=( $(compgen -W "> a11,a12,a13,a21,a22,a23,a31,a32,a33") ); return 0;;
        "mmul" | "-mmul" | "+mmul")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "mod" | "-mod" | "+mod")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "mode3d" | "-mode3d" | "+mode3d")
            COMPREPLY=( $(compgen -W "> _mode") ); return 0;;
        "moded3d" | "-moded3d" | "+moded3d")
            COMPREPLY=( $(compgen -W "> _mode") ); return 0;;
        "montage" | "-montage" | "+montage")
            COMPREPLY=( $(compgen -W "> \"_layout_code\",_montage_mode={_0<=centering<=1_|_2<=scale+2<=3_},_output_mode={_0=single_layer_|_1=multiple_layers_},\"_processing_command\"") ); return 0;;
        "morph" | "-morph" | "+morph")
            COMPREPLY=( $(compgen -W "> nb_inner_frames>=1,_smoothness>=0,_precision>=0") ); return 0;;
        "morph_files" | "-morph_files" | "+morph_files")
            COMPREPLY=( $(compgen -W "> \"filename_pattern\",_nb_inner_frames>0,_smoothness>=0,_precision>=0,_first_frame>=0,_last_frame={_>=0_|_-1=last_},_frame_step>=1,_output_filename") ); return 0;;
        "morph_rbf" | "-morph_rbf" | "+morph_rbf")
            COMPREPLY=( $(compgen -W "> nb_inner_frames>=1,xs0[%],ys0[%],xt0[%],yt0[%],...,xsN[%],ysN[%],xtN[%],ytN[%]") ); return 0;;
        "morph_video" | "-morph_video" | "+morph_video")
            COMPREPLY=( $(compgen -W "> video_filename,_nb_inner_frames>0,_smoothness>=0,_precision>=0,_first_frame>=0,_last_frame={_>=0_|_-1=last_},_frame_step>=1,_output_filename") ); return 0;;
        "mosaic" | "-mosaic" | "+mosaic")
            COMPREPLY=( $(compgen -W "> 0<=_density<=100") ); return 0;;
        "move" | "-move" | "+move")
            COMPREPLY=( $(compgen -W "> position[%]") ); return 0;;
        "mproj" | "-mproj" | "+mproj")
            COMPREPLY=( $(compgen -W "> [dictionary],_method,_max_iter={_0=auto_|_>0_},_max_residual>=0") ); return 0;;
        "mse" | "-mse" | "+mse")
            COMPREPLY=( $(compgen -W "> [reference]") ); return 0;;
        "mul3d" | "-mul3d" | "+mul3d")
            COMPREPLY=( $(compgen -W "factor factor_x,factor_y,_factor_z") ); return 0;;
        "mul" | "-mul" | "+mul")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "mul_channels" | "-mul_channels" | "+mul_channels")
            COMPREPLY=( $(compgen -W "> value1,_value2,...,_valueN") ); return 0;;
        "mul_complex" | "-mul_complex" | "+mul_complex")
            COMPREPLY=( $(compgen -W "> [multiplier_real,multiplier_imag]") ); return 0;;
        "mutex" | "-mutex" | "+mutex")
            COMPREPLY=( $(compgen -W "> index,_action={_0=unlock_|_1=lock_}") ); return 0;;
        "mv" | "-mv" | "+mv")
            COMPREPLY=( $(compgen -W "> position[%]") ); return 0;;
        "n" | "-n" | "+n")
            COMPREPLY=( $(compgen -W "{_value0[%]_|_[image0]_},{_value1[%]_|_[image1]_},_constant_case_ratio [image]") ); return 0;;
        "name" | "-name" | "+name")
            COMPREPLY=( $(compgen -W "> \"name1\",\"name2\",...") ); return 0;;
        "named" | "-named" | "+named")
            COMPREPLY=( $(compgen -W "> _mode,\"name1\",\"name2\",...") ); return 0;;
        "negate" | "-negate" | "+negate")
            COMPREPLY=( $(compgen -W "base_value (no_arg)") ); return 0;;
        "neq" | "-neq" | "+neq")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "network" | "-network" | "+network")
            COMPREPLY=( $(compgen -W "> mode={_-1=disabled_|_0=enabled_w/o_timeout_|_>0=enabled_w/_specified_timeout_in_seconds_}") ); return 0;;
        "newton_fractal" | "-newton_fractal" | "+newton_fractal")
            COMPREPLY=( $(compgen -W "> z0r,z0i,z1r,z1i,_angle,0<=_descent_method<=2,_iteration_max>=0,_convergence_precision>0,_expr_p(z),_expr_dp(z),_expr_d2p(z)") ); return 0;;
        "nlmeans" | "-nlmeans" | "+nlmeans")
            COMPREPLY=( $(compgen -W "[guide],_patch_radius>0,_spatial_bandwidth>0,_tonal_bandwidth>0,_patch_measure_command _patch_radius>0,_spatial_bandwidth>0,_tonal_bandwidth>0,_patch_measure_command") ); return 0;;
        "nlmeans_core" | "-nlmeans_core" | "+nlmeans_core")
            COMPREPLY=( $(compgen -W "> _reference_image,_scaling_map,_patch_radius>0,_spatial_bandwidth>0") ); return 0;;
        "nm" | "-nm" | "+nm")
            COMPREPLY=( $(compgen -W "> \"name1\",\"name2\",...") ); return 0;;
        "nmd" | "-nmd" | "+nmd")
            COMPREPLY=( $(compgen -W "> _mode,\"name1\",\"name2\",...") ); return 0;;
        "nn_check_layer" | "-nn_check_layer" | "+nn_check_layer")
            COMPREPLY=( $(compgen -W "> name") ); return 0;;
        "nn_layer_add" | "-nn_layer_add" | "+nn_layer_add")
            COMPREPLY=( $(compgen -W "> name,in0,in1") ); return 0;;
        "nn_layer_append" | "-nn_layer_append" | "+nn_layer_append")
            COMPREPLY=( $(compgen -W "> name,in0,in1") ); return 0;;
        "nn_layer_avgpool2d" | "-nn_layer_avgpool2d" | "+nn_layer_avgpool2d")
            COMPREPLY=( $(compgen -W "> name,in") ); return 0;;
        "nn_layer_batchnorm" | "-nn_layer_batchnorm" | "+nn_layer_batchnorm")
            COMPREPLY=( $(compgen -W "> name,in,_learning_mode.") ); return 0;;
        "nn_layer_clone" | "-nn_layer_clone" | "+nn_layer_clone")
            COMPREPLY=( $(compgen -W "> name0,name1,in") ); return 0;;
        "nn_layer_conv2d" | "-nn_layer_conv2d" | "+nn_layer_conv2d")
            COMPREPLY=( $(compgen -W "> name,in,nb_channels>0,_kernel_size>0,_stride>0,_dilation,_is_learned={_0_|_1_}") ); return 0;;
        "nn_layer_conv2dbnnl" | "-nn_layer_conv2dbnnl" | "+nn_layer_conv2dbnnl")
            COMPREPLY=( $(compgen -W "> name,in,nb_channels>0,_kernel_size>0,_stride>0,_dilation>0,_activation,_is_learned={_0_|_1_}") ); return 0;;
        "nn_layer_conv2dnl" | "-nn_layer_conv2dnl" | "+nn_layer_conv2dnl")
            COMPREPLY=( $(compgen -W "> name,in,nb_channels>0,_kernel_size>0,_stride>0,_dilation>0,_activation,_is_learned={_0_|_1_}") ); return 0;;
        "nn_layer_crop" | "-nn_layer_crop" | "+nn_layer_crop")
            COMPREPLY=( $(compgen -W "> name,in,x0,y0,z0,x1,y1,z1") ); return 0;;
        "nn_layer_fc" | "-nn_layer_fc" | "+nn_layer_fc")
            COMPREPLY=( $(compgen -W "> name,in,nb_channels>0,_is_learned={_0_|_1_}") ); return 0;;
        "nn_layer_fcbnnl" | "-nn_layer_fcbnnl" | "+nn_layer_fcbnnl")
            COMPREPLY=( $(compgen -W "> name,in,nb_neurons>0,_activation,_is_learned={_0_|_1_}") ); return 0;;
        "nn_layer_fcnl" | "-nn_layer_fcnl" | "+nn_layer_fcnl")
            COMPREPLY=( $(compgen -W "> name,in,nb_neurons>0,_activation,_is_learned={_0_|_1_}") ); return 0;;
        "nn_layer_input" | "-nn_layer_input" | "+nn_layer_input")
            COMPREPLY=( $(compgen -W "> name,width,_height,_depth,_spectrum") ); return 0;;
        "nn_layer_maxpool2d" | "-nn_layer_maxpool2d" | "+nn_layer_maxpool2d")
            COMPREPLY=( $(compgen -W "> name,in") ); return 0;;
        "nn_layer_nl" | "-nn_layer_nl" | "+nn_layer_nl")
            COMPREPLY=( $(compgen -W "> name,in,_activation") ); return 0;;
        "nn_layer_rename" | "-nn_layer_rename" | "+nn_layer_rename")
            COMPREPLY=( $(compgen -W "> name,in") ); return 0;;
        "nn_layer_resconv2d" | "-nn_layer_resconv2d" | "+nn_layer_resconv2d")
            COMPREPLY=( $(compgen -W "> name,in,_kernel_size>0,_stride>0,_dilation>0,_is_learned={_0_|_1_}") ); return 0;;
        "nn_layer_resconv2dnl" | "-nn_layer_resconv2dnl" | "+nn_layer_resconv2dnl")
            COMPREPLY=( $(compgen -W "> name,in,_kernel_size>0,_stride>0,_dilation>0,_activation,_is_learned={_0_|_1_}") ); return 0;;
        "nn_layer_reshape" | "-nn_layer_reshape" | "+nn_layer_reshape")
            COMPREPLY=( $(compgen -W "> name,in,width>0,height>0,depth>0,spectrum>0") ); return 0;;
        "nn_layer_resize" | "-nn_layer_resize" | "+nn_layer_resize")
            COMPREPLY=( $(compgen -W "> name,in,width[%]>0,_height[%]>0,_depth[%]>0,_spectrum[%]>0,_interpolation") ); return 0;;
        "nn_layer_run" | "-nn_layer_run" | "+nn_layer_run")
            COMPREPLY=( $(compgen -W "> name,in,\"command\",_width[%]>0,_height[%]>0,_depth[%]>0,_spectrum[%]>0") ); return 0;;
        "nn_layer_split" | "-nn_layer_split" | "+nn_layer_split")
            COMPREPLY=( $(compgen -W "> name0,name1,in,nb_channels0") ); return 0;;
        "nn_load" | "-nn_load" | "+nn_load")
            COMPREPLY=( $(compgen -W "> \'filename.gmz\'") ); return 0;;
        "nn_loss_bce" | "-nn_loss_bce" | "+nn_loss_bce")
            COMPREPLY=( $(compgen -W "> name,in,ground_truth") ); return 0;;
        "nn_loss_mse" | "-nn_loss_mse" | "+nn_loss_mse")
            COMPREPLY=( $(compgen -W "> name,in,ground_truth") ); return 0;;
        "nn_save" | "-nn_save" | "+nn_save")
            COMPREPLY=( $(compgen -W "> \'filename.gmz\'") ); return 0;;
        "nn_trainer" | "-nn_trainer" | "+nn_trainer")
            COMPREPLY=( $(compgen -W "> name,loss,_learning_rate>0,_optimizer,_scheduler") ); return 0;;
        "noise" | "-noise" | "+noise")
            COMPREPLY=( $(compgen -W "> std_deviation>=0[%],_noise_type") ); return 0;;
        "noise_hurl" | "-noise_hurl" | "+noise_hurl")
            COMPREPLY=( $(compgen -W "> _amplitude>=0") ); return 0;;
        "noise_perlin" | "-noise_perlin" | "+noise_perlin")
            COMPREPLY=( $(compgen -W "> _scale_x[%]>0,_scale_y[%]>0,_scale_z[%]>0,_seed_x,_seed_y,_seed_z") ); return 0;;
        "noise_poissondisk" | "-noise_poissondisk" | "+noise_poissondisk")
            COMPREPLY=( $(compgen -W "> _radius[%]>0,_max_sample_attempts>0,_p_norm>0") ); return 0;;
        "normalize" | "-normalize" | "+normalize")
            COMPREPLY=( $(compgen -W "{_value0[%]_|_[image0]_},{_value1[%]_|_[image1]_},_constant_case_ratio [image]") ); return 0;;
        "normalize_filename" | "-normalize_filename" | "+normalize_filename")
            COMPREPLY=( $(compgen -W "> filename") ); return 0;;
        "normalize_local" | "-normalize_local" | "+normalize_local")
            COMPREPLY=( $(compgen -W "> _amplitude>=0,_radius>0,_n_smooth>=0[%],_a_smooth>=0[%],_is_cut={_0_|_1_},_min=0,_max=255") ); return 0;;
        "normalized_cross_correlation" | "-normalized_cross_correlation" | "+normalized_cross_correlation")
            COMPREPLY=( $(compgen -W "> [mask]") ); return 0;;
        "normp" | "-normp" | "+normp")
            COMPREPLY=( $(compgen -W "> p>=0") ); return 0;;
        "o3d" | "-o3d" | "+o3d")
            COMPREPLY=( $(compgen -W "> _opacity") ); return 0;;
        "object3d" | "-object3d" | "+object3d")
            COMPREPLY=( $(compgen -W "> [object3d],_x[%],_y[%],_z,_opacity,_rendering_mode,_is_double_sided={_0_|_1_},_is_zbuffer={_0_|_1_},_focale,_light_x,_light_y,_light_z,_specular_lightness,_specular_shininess") ); return 0;;
        "oct2dec" | "-oct2dec" | "+oct2dec")
            COMPREPLY=( $(compgen -W "> octal_int1,...") ); return 0;;
        "oct" | "-oct" | "+oct")
            COMPREPLY=( $(compgen -W "> octal_int1,...") ); return 0;;
        "on" | "-on" | "+on")
            COMPREPLY=( $(compgen -W "> filename,_index") ); return 0;;
        "op" | "-op" | "+op")
            COMPREPLY=( $(compgen -W "> prefix") ); return 0;;
        "opacity3d" | "-opacity3d" | "+opacity3d")
            COMPREPLY=( $(compgen -W "> _opacity") ); return 0;;
        "opening" | "-opening" | "+opening")
            COMPREPLY=( $(compgen -W "size>=0 size_x>=0,size_y>=0,_size_z>=0 [kernel],_boundary_conditions,_is_real={_0=binary-mode_|_1=real-mode_}") ); return 0;;
        "opening_circ" | "-opening_circ" | "+opening_circ")
            COMPREPLY=( $(compgen -W "> _size>=0,_is_real={_0_|_1_}") ); return 0;;
        "or" | "-or" | "+or")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "orthogonalize" | "-orthogonalize" | "+orthogonalize")
            COMPREPLY=( $(compgen -W "> _mode_=_{_0=orthogonalize_|_1=orthonormalize_}") ); return 0;;
        "otsu" | "-otsu" | "+otsu")
            COMPREPLY=( $(compgen -W "> _nb_levels>0") ); return 0;;
        "ox" | "-ox" | "+ox")
            COMPREPLY=( $(compgen -W "> extension1,_extension2,_...,_extensionN,_output_at_same_location={_0_|_1_}") ); return 0;;
        "p3d" | "-p3d" | "+p3d")
            COMPREPLY=( $(compgen -W "> mode") ); return 0;;
        "pack" | "-pack" | "+pack")
            COMPREPLY=( $(compgen -W "> is_ratio_constraint={_0_|_1_},_sort_criterion") ); return 0;;
        "pack_sprites" | "-pack_sprites" | "+pack_sprites")
            COMPREPLY=( $(compgen -W "> _nb_scales>=0,0<=_min_scale<=100,_allow_rotation={_0=0_deg._|_1=180_deg._|_2=90_deg._|_3=any_},_spacing,_precision>=0,max_iterations>=0") ); return 0;;
        "padint" | "-padint" | "+padint")
            COMPREPLY=( $(compgen -W "> number,_size>0") ); return 0;;
        "palette" | "-palette" | "+palette")
            COMPREPLY=( $(compgen -W "> palette_name_|_palette_number") ); return 0;;
        "parallel" | "-parallel" | "+parallel")
            COMPREPLY=( $(compgen -W "> _wait_threads,\"command1\",\"command2\",...") ); return 0;;
        "parametric3d" | "-parametric3d" | "+parametric3d")
            COMPREPLY=( $(compgen -W "> _x(a,b),_y(a,b),_z(a,b),_amin,_amax,_bmin,_bmax,_res_a>0,_res_b>0,_res_x>0,_res_y>0,_res_z>0,_smoothness>=0,_isovalue>=0") ); return 0;;
        "parse_cli" | "-parse_cli" | "+parse_cli")
            COMPREPLY=( $(compgen -W "> _output_mode,_{_*_|_command_name_}") ); return 0;;
        "parse_gui" | "-parse_gui" | "+parse_gui")
            COMPREPLY=( $(compgen -W "> _outputmode,_{_*_|_filter_name}") ); return 0;;
        "pass" | "-pass" | "+pass")
            COMPREPLY=( $(compgen -W "> _shared_state={_-1=status_only_|_0=non-shared_(copy)_|_1=shared_|_2=adaptive_}") ); return 0;;
        "patches2img" | "-patches2img" | "+patches2img")
            COMPREPLY=( $(compgen -W "> width>0,height>0,_overlap[%]>0,_overlap_std[%]") ); return 0;;
        "patches" | "-patches" | "+patches")
            COMPREPLY=( $(compgen -W "> patch_width>0,patch_height>0,patch_depth>0,x0,y0,z0,_x1,_y1,_z1,...,_xN,_yN,_zN") ); return 0;;
        "pca_patch3d" | "-pca_patch3d" | "+pca_patch3d")
            COMPREPLY=( $(compgen -W "> _patch_size>0,_M>0,_N>0,_normalize_input={_0_|_1_},_normalize_output={_0_|_1_},_lambda_xy") ); return 0;;
        "pde_flow" | "-pde_flow" | "+pde_flow")
            COMPREPLY=( $(compgen -W "> _nb_iter>=0,_dt,_velocity_command,_keep_sequence={_0_|_1_}") ); return 0;;
        "pencilbw" | "-pencilbw" | "+pencilbw")
            COMPREPLY=( $(compgen -W "> _size>=0,_amplitude>=0") ); return 0;;
        "percentile" | "-percentile" | "+percentile")
            COMPREPLY=( $(compgen -W "> [mask],0<=_min_percentile[%]<=100,0<=_max_percentile[%]<=100.") ); return 0;;
        "permute" | "-permute" | "+permute")
            COMPREPLY=( $(compgen -W "> permutation_string") ); return 0;;
        "peronamalik_flow" | "-peronamalik_flow" | "+peronamalik_flow")
            COMPREPLY=( $(compgen -W "> K_factor>0,_nb_iter>=0,_dt,_keep_sequence={_0_|_1_}") ); return 0;;
        "phase_correlation" | "-phase_correlation" | "+phase_correlation")
            COMPREPLY=( $(compgen -W "> [destination]") ); return 0;;
        "piechart" | "-piechart" | "+piechart")
            COMPREPLY=( $(compgen -W "> label_height>=0,label_R,label_G,label_B,\"label1\",value1,R1,G1,B1,...,\"labelN\",valueN,RN,GN,BN") ); return 0;;
        "pixelize" | "-pixelize" | "+pixelize")
            COMPREPLY=( $(compgen -W "> _scale_x>0,_scale_y>0,_scale_z>0") ); return 0;;
        "pixelsort" | "-pixelsort" | "+pixelsort")
            COMPREPLY=( $(compgen -W "> _ordering={_+_|_-_},_axis={_x_|_y_|_z_|_xy_|_yx_},_[sorting_criterion],_[mask]") ); return 0;;
        "plane3d" | "-plane3d" | "+plane3d")
            COMPREPLY=( $(compgen -W "> _size_x,_size_y,_nb_subdivisions_x>0,_nb_subdisivions_y>0") ); return 0;;
        "plasma" | "-plasma" | "+plasma")
            COMPREPLY=( $(compgen -W "> _alpha,_beta,_scale>=0") ); return 0;;
        "plot" | "-plot" | "+plot")
            COMPREPLY=( $(compgen -W "_plot_type,_vertex_type,_xmin,_xmax,_ymin,_ymax,_exit_on_anykey={_0_|_1_} \'formula\',_resolution>=0,_plot_type,_vertex_type,_xmin,xmax,_ymin,_ymax,_exit_on_anykey={_0_|_1_}") ); return 0;;
        "point3d" | "-point3d" | "+point3d")
            COMPREPLY=( $(compgen -W "> x0,y0,z0") ); return 0;;
        "point" | "-point" | "+point")
            COMPREPLY=( $(compgen -W "> x[%],_y[%],_z[%],_opacity,_color1,...") ); return 0;;
        "pointcloud" | "-pointcloud" | "+pointcloud")
            COMPREPLY=( $(compgen -W "> _type_=_{_-X=-X-opacity_|_0=binary_|_1=cumulative_|_2=label_|_3=retrieve_coordinates_},_width,_height>0,_depth>0") ); return 0;;
        "polar2euclidean" | "-polar2euclidean" | "+polar2euclidean")
            COMPREPLY=( $(compgen -W "> _center_x[%],_center_y[%],_stretch_factor>0,_boundary_conditions={_0=dirichlet_|_1=neumann_|_2=periodic_|_3=mirror_}") ); return 0;;
        "polaroid" | "-polaroid" | "+polaroid")
            COMPREPLY=( $(compgen -W "> _size1>=0,_size2>=0") ); return 0;;
        "polka_dots" | "-polka_dots" | "+polka_dots")
            COMPREPLY=( $(compgen -W "> diameter>=0,_density,_offset1,_offset2,_angle,_aliasing,_shading,_opacity,_color,...") ); return 0;;
        "polygon" | "-polygon" | "+polygon")
            COMPREPLY=( $(compgen -W "> N>=1,x1[%],y1[%],...,xN[%],yN[%],_opacity,_pattern,_color1,...") ); return 0;;
        "polygonize" | "-polygonize" | "+polygonize")
            COMPREPLY=( $(compgen -W "> _warp_amplitude>=0,_smoothness[%]>=0,_min_area[%]>=0,_resolution_x[%]>0,_resolution_y[%]>0") ); return 0;;
        "portrait" | "-portrait" | "+portrait")
            COMPREPLY=( $(compgen -W "> _size>0") ); return 0;;
        "pose3d" | "-pose3d" | "+pose3d")
            COMPREPLY=( $(compgen -W "> p1,...,p12") ); return 0;;
        "poster_edges" | "-poster_edges" | "+poster_edges")
            COMPREPLY=( $(compgen -W "> 0<=_edge_threshold<=100,0<=_edge_shade<=100,_edge_thickness>=0,_edge_antialiasing>=0,0<=_posterization_level<=15,_posterization_antialiasing>=0") ); return 0;;
        "poster_hope" | "-poster_hope" | "+poster_hope")
            COMPREPLY=( $(compgen -W "> _smoothness>=0") ); return 0;;
        "pow" | "-pow" | "+pow")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "primitives3d" | "-primitives3d" | "+primitives3d")
            COMPREPLY=( $(compgen -W "> mode") ); return 0;;
        "progress" | "-progress" | "+progress")
            COMPREPLY=( $(compgen -W "0<=value<=100 -1") ); return 0;;
        "projections3d" | "-projections3d" | "+projections3d")
            COMPREPLY=( $(compgen -W "> _x[%],_y[%],_z[%],_is_bounding_box={_0_|_1_}") ); return 0;;
        "pseudogray" | "-pseudogray" | "+pseudogray")
            COMPREPLY=( $(compgen -W "> _max_increment>=0,_JND_threshold>=0,_bits_depth>0") ); return 0;;
        "psnr" | "-psnr" | "+psnr")
            COMPREPLY=( $(compgen -W "> [reference],_max_value>0") ); return 0;;
        "psnr_matrix" | "-psnr_matrix" | "+psnr_matrix")
            COMPREPLY=( $(compgen -W "> _max_value>0") ); return 0;;
        "puzzle" | "-puzzle" | "+puzzle")
            COMPREPLY=( $(compgen -W "> _width>0,_height>0,_M>=1,_N>=1,_curvature,_centering,_connectors_variability,_resolution>=1") ); return 0;;
        "pyramid3d" | "-pyramid3d" | "+pyramid3d")
            COMPREPLY=( $(compgen -W "> width,height") ); return 0;;
        "quadrangle3d" | "-quadrangle3d" | "+quadrangle3d")
            COMPREPLY=( $(compgen -W "> x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3") ); return 0;;
        "quadratize_tiles" | "-quadratize_tiles" | "+quadratize_tiles")
            COMPREPLY=( $(compgen -W "> M>0,_N>0") ); return 0;;
        "quantize" | "-quantize" | "+quantize")
            COMPREPLY=( $(compgen -W "> nb_levels>=1,_keep_values={_0_|_1_},_quantization_type={_-1=median-cut_|_0=k-means_|_1=uniform_}") ); return 0;;
        "quantize_area" | "-quantize_area" | "+quantize_area")
            COMPREPLY=( $(compgen -W "> _min_area>0") ); return 0;;
        "quiver" | "-quiver" | "+quiver")
            COMPREPLY=( $(compgen -W "> [function_image],_sampling[%]>0,_factor>=0,_is_arrow={_0_|_1_},_opacity,_color1,...") ); return 0;;
        "r2din" | "-r2din" | "+r2din")
            COMPREPLY=( $(compgen -W "> width[%]>0,_height[%]>0,_interpolation,_boundary_conditions,_ax,_ay,_az,_ac") ); return 0;;
        "r2dout" | "-r2dout" | "+r2dout")
            COMPREPLY=( $(compgen -W "> width[%]>0,_height[%]>0,_interpolation,_boundary_conditions,_ax,_ay,_az,_ac") ); return 0;;
        "r2dx" | "-r2dx" | "+r2dx")
            COMPREPLY=( $(compgen -W "> width[%]>0,_interpolation,_boundary_conditions,_ax,_ay,_az,_ac") ); return 0;;
        "r2dy" | "-r2dy" | "+r2dy")
            COMPREPLY=( $(compgen -W "> height[%]>=0,_interpolation,_boundary_conditions,_ax,_ay,_az,_ac") ); return 0;;
        "r3d" | "-r3d" | "+r3d")
            COMPREPLY=( $(compgen -W "> u,v,w,angle") ); return 0;;
        "r3din" | "-r3din" | "+r3din")
            COMPREPLY=( $(compgen -W "> width[%]>0,_height[%]>0,_depth[%]>0,_interpolation,_boundary_conditions,_ax,_ay,_az,_ac") ); return 0;;
        "r3dout" | "-r3dout" | "+r3dout")
            COMPREPLY=( $(compgen -W "> width[%]>0,_height[%]>0,_depth[%]>0,_interpolation,_boundary_conditions,_ax,_ay,_az,_ac") ); return 0;;
        "r3dx" | "-r3dx" | "+r3dx")
            COMPREPLY=( $(compgen -W "> width[%]>0,_interpolation,_boundary_conditions,_ax,_ay,_az,_ac") ); return 0;;
        "r3dy" | "-r3dy" | "+r3dy")
            COMPREPLY=( $(compgen -W "> height[%]>0,_interpolation,_boundary_conditions,_ax,_ay,_az,_ac") ); return 0;;
        "r3dz" | "-r3dz" | "+r3dz")
            COMPREPLY=( $(compgen -W "> depth[%]>0,_interpolation,_boundary_conditions,_ax,_ay,_az,_ac") ); return 0;;
        "r" | "-r" | "+r")
            COMPREPLY=( $(compgen -W "> {[image_w]_|_width>0[%]},_{[image_h]_|_height>0[%]},_{[image_d]_|_depth>0[%]},_{[image_s]_|_spectrum>0[%]},_interpolation,_boundary_conditions,_ax,_ay,_az,_ac") ); return 0;;
        "raindrops" | "-raindrops" | "+raindrops")
            COMPREPLY=( $(compgen -W "> _amplitude,_density>=0,_wavelength>=0,_merging_steps>=0") ); return 0;;
        "rand" | "-rand" | "+rand")
            COMPREPLY=( $(compgen -W "{_value0[%]_|_[image0]_},_{_value1[%]_|_[image1]_} [image]") ); return 0;;
        "random3d" | "-random3d" | "+random3d")
            COMPREPLY=( $(compgen -W "> nb_points>=0") ); return 0;;
        "random_pattern" | "-random_pattern" | "+random_pattern")
            COMPREPLY=( $(compgen -W "> _width>0,_height>0,_min_detail_level>=0") ); return 0;;
        "rbf" | "-rbf" | "+rbf")
            COMPREPLY=( $(compgen -W "dx,_x0,_x1,_phi(r) dx,dy,_x0,_y0,_x1,_y1,_phi(r) dx,dy,dz,x0,y0,z0,x1,y1,z1,phi(r)") ); return 0;;
        "rectangle" | "-rectangle" | "+rectangle")
            COMPREPLY=( $(compgen -W "> x0[%],y0[%],x1[%],y1[%],_opacity,_pattern,_color1,...") ); return 0;;
        "red_eye" | "-red_eye" | "+red_eye")
            COMPREPLY=( $(compgen -W "> 0<=_threshold<=100,_smoothness>=0,0<=attenuation<=1") ); return 0;;
        "register_nonrigid" | "-register_nonrigid" | "+register_nonrigid")
            COMPREPLY=( $(compgen -W "> [destination],_smoothness>=0,_precision>0,_nb_scale>=0") ); return 0;;
        "register_rigid" | "-register_rigid" | "+register_rigid")
            COMPREPLY=( $(compgen -W "> [destination],_smoothness>=0,_boundary_conditions={_0=dirichlet_|_1=neumann_|_2=periodic_|_3=mirror_}") ); return 0;;
        "remove_copymark" | "-remove_copymark" | "+remove_copymark")
            COMPREPLY=( $(compgen -W "> \"image_name\"") ); return 0;;
        "remove_hotpixels" | "-remove_hotpixels" | "+remove_hotpixels")
            COMPREPLY=( $(compgen -W "> _mask_size>0,__threshold[%]>0") ); return 0;;
        "remove_named" | "-remove_named" | "+remove_named")
            COMPREPLY=( $(compgen -W "> \"name1\",\"name2\",...") ); return 0;;
        "remove_pixels" | "-remove_pixels" | "+remove_pixels")
            COMPREPLY=( $(compgen -W "> number_of_pixels[%]>=0") ); return 0;;
        "repeat" | "-repeat" | "+repeat")
            COMPREPLY=( $(compgen -W "> nb_iterations") ); return 0;;
        "replace" | "-replace" | "+replace")
            COMPREPLY=( $(compgen -W "> source,target") ); return 0;;
        "replace_color" | "-replace_color" | "+replace_color")
            COMPREPLY=( $(compgen -W "> tolerance[%]>=0,smoothness[%]>=0,src1,src2,...,dest1,dest2,...") ); return 0;;
        "replace_inf" | "-replace_inf" | "+replace_inf")
            COMPREPLY=( $(compgen -W "> _expression") ); return 0;;
        "replace_nan" | "-replace_nan" | "+replace_nan")
            COMPREPLY=( $(compgen -W "> _expression") ); return 0;;
        "replace_naninf" | "-replace_naninf" | "+replace_naninf")
            COMPREPLY=( $(compgen -W "> _expression") ); return 0;;
        "replace_seq" | "-replace_seq" | "+replace_seq")
            COMPREPLY=( $(compgen -W "> \"search_seq\",\"replace_seq\"") ); return 0;;
        "replace_str" | "-replace_str" | "+replace_str")
            COMPREPLY=( $(compgen -W "> \"search_str\",\"replace_str\"") ); return 0;;
        "resize2din" | "-resize2din" | "+resize2din")
            COMPREPLY=( $(compgen -W "> width[%]>0,_height[%]>0,_interpolation,_boundary_conditions,_ax,_ay,_az,_ac") ); return 0;;
        "resize2dout" | "-resize2dout" | "+resize2dout")
            COMPREPLY=( $(compgen -W "> width[%]>0,_height[%]>0,_interpolation,_boundary_conditions,_ax,_ay,_az,_ac") ); return 0;;
        "resize2dx" | "-resize2dx" | "+resize2dx")
            COMPREPLY=( $(compgen -W "> width[%]>0,_interpolation,_boundary_conditions,_ax,_ay,_az,_ac") ); return 0;;
        "resize2dy" | "-resize2dy" | "+resize2dy")
            COMPREPLY=( $(compgen -W "> height[%]>=0,_interpolation,_boundary_conditions,_ax,_ay,_az,_ac") ); return 0;;
        "resize3din" | "-resize3din" | "+resize3din")
            COMPREPLY=( $(compgen -W "> width[%]>0,_height[%]>0,_depth[%]>0,_interpolation,_boundary_conditions,_ax,_ay,_az,_ac") ); return 0;;
        "resize3dout" | "-resize3dout" | "+resize3dout")
            COMPREPLY=( $(compgen -W "> width[%]>0,_height[%]>0,_depth[%]>0,_interpolation,_boundary_conditions,_ax,_ay,_az,_ac") ); return 0;;
        "resize3dx" | "-resize3dx" | "+resize3dx")
            COMPREPLY=( $(compgen -W "> width[%]>0,_interpolation,_boundary_conditions,_ax,_ay,_az,_ac") ); return 0;;
        "resize3dy" | "-resize3dy" | "+resize3dy")
            COMPREPLY=( $(compgen -W "> height[%]>0,_interpolation,_boundary_conditions,_ax,_ay,_az,_ac") ); return 0;;
        "resize3dz" | "-resize3dz" | "+resize3dz")
            COMPREPLY=( $(compgen -W "> depth[%]>0,_interpolation,_boundary_conditions,_ax,_ay,_az,_ac") ); return 0;;
        "resize" | "-resize" | "+resize")
            COMPREPLY=( $(compgen -W "> {[image_w]_|_width>0[%]},_{[image_h]_|_height>0[%]},_{[image_d]_|_depth>0[%]},_{[image_s]_|_spectrum>0[%]},_interpolation,_boundary_conditions,_ax,_ay,_az,_ac") ); return 0;;
        "resize_as_image" | "-resize_as_image" | "+resize_as_image")
            COMPREPLY=( $(compgen -W "> [reference],_interpolation,_boundary_conditions,_ax,_ay,_az,_ac") ); return 0;;
        "resize_mn" | "-resize_mn" | "+resize_mn")
            COMPREPLY=( $(compgen -W "> width[%]>=0,_height[%]>=0,_depth[%]>=0,_B_value,_C_value") ); return 0;;
        "resize_pow2" | "-resize_pow2" | "+resize_pow2")
            COMPREPLY=( $(compgen -W "> _interpolation,_boundary_conditions,_ax,_ay,_az,_ac") ); return 0;;
        "resize_ratio2d" | "-resize_ratio2d" | "+resize_ratio2d")
            COMPREPLY=( $(compgen -W "> width>0,height>0,_mode={_0=inside_|_1=outside_|_2=padded_},0=<_interpolation<=6") ); return 0;;
        "retinex" | "-retinex" | "+retinex")
            COMPREPLY=( $(compgen -W "> _value_offset>0,_colorspace={_hsi_|_hsv_|_lab_|_lrgb_|_rgb_|_ycbcr_},0<=_min_cut<=100,0<=_max_cut<=100,_sigma_low>0,_sigma_mid>0,_sigma_high>0") ); return 0;;
        "rgb2bayer" | "-rgb2bayer" | "+rgb2bayer")
            COMPREPLY=( $(compgen -W "> _start_pattern=0,_color_grid=0") ); return 0;;
        "rgb2jzazbz" | "-rgb2jzazbz" | "+rgb2jzazbz")
            COMPREPLY=( $(compgen -W "illuminant={_0=D50_|_1=D65_|_2=E_} (no_arg)") ); return 0;;
        "rgb2lab8" | "-rgb2lab8" | "+rgb2lab8")
            COMPREPLY=( $(compgen -W "illuminant={_0=D50_|_1=D65_|_2=E_} (no_arg)") ); return 0;;
        "rgb2lab" | "-rgb2lab" | "+rgb2lab")
            COMPREPLY=( $(compgen -W "illuminant={_0=D50_|_1=D65_|_2=E_} (no_arg)") ); return 0;;
        "rgb2lch8" | "-rgb2lch8" | "+rgb2lch8")
            COMPREPLY=( $(compgen -W "illuminant={_0=D50_|_1=D65_|_2=E_} (no_arg)") ); return 0;;
        "rgb2lch" | "-rgb2lch" | "+rgb2lch")
            COMPREPLY=( $(compgen -W "illuminant={_0=D50_|_1=D65_|_2=E_} (no_arg)") ); return 0;;
        "rgb2xyz8" | "-rgb2xyz8" | "+rgb2xyz8")
            COMPREPLY=( $(compgen -W "illuminant={_0=D50_|_1=D65_|_2=E_} (no_arg)") ); return 0;;
        "rgb2xyz" | "-rgb2xyz" | "+rgb2xyz")
            COMPREPLY=( $(compgen -W "illuminant={_0=D50_|_1=D65_|_2=E_} (no_arg)") ); return 0;;
        "ri" | "-ri" | "+ri")
            COMPREPLY=( $(compgen -W "> [reference],_interpolation,_boundary_conditions,_ax,_ay,_az,_ac") ); return 0;;
        "ripple" | "-ripple" | "+ripple")
            COMPREPLY=( $(compgen -W "> _amplitude,_bandwidth,_shape={_0=block_|_1=triangle_|_2=sine_|_3=sine+_|_4=random_},_angle,_offset") ); return 0;;
        "rmn" | "-rmn" | "+rmn")
            COMPREPLY=( $(compgen -W "> \"name1\",\"name2\",...") ); return 0;;
        "rodilius" | "-rodilius" | "+rodilius")
            COMPREPLY=( $(compgen -W "> 0<=_amplitude<=100,_0<=thickness<=100,_sharpness>=0,_nb_orientations>0,_offset,_color_mode={_0=darker_|_1=brighter_}") ); return 0;;
        "rol" | "-rol" | "+rol")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "rolling_guidance" | "-rolling_guidance" | "+rolling_guidance")
            COMPREPLY=( $(compgen -W "> std_deviation_s[%]>=0,std_deviation_r[%]>=0,_precision>=0") ); return 0;;
        "ror" | "-ror" | "+ror")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "rorschach" | "-rorschach" | "+rorschach")
            COMPREPLY=( $(compgen -W "> \'smoothness[%]>=0\',\'mirroring={_0=none_|_1=x_|_2=y_|_3=xy_}") ); return 0;;
        "rotate3d" | "-rotate3d" | "+rotate3d")
            COMPREPLY=( $(compgen -W "> u,v,w,angle") ); return 0;;
        "rotate" | "-rotate" | "+rotate")
            COMPREPLY=( $(compgen -W "angle,_interpolation,_boundary_conditions,_center_x[%],_center_y[%] u,v,w,angle,interpolation,boundary_conditions,_center_x[%],_center_y[%],_center_z[%]") ); return 0;;
        "rotate_tileable" | "-rotate_tileable" | "+rotate_tileable")
            COMPREPLY=( $(compgen -W "> angle,_max_size_factor>=0") ); return 0;;
        "rotate_tiles" | "-rotate_tiles" | "+rotate_tiles")
            COMPREPLY=( $(compgen -W "> angle,_M>0,N>0") ); return 0;;
        "rotation3d" | "-rotation3d" | "+rotation3d")
            COMPREPLY=( $(compgen -W "> u,v,w,angle") ); return 0;;
        "rotoidoscope" | "-rotoidoscope" | "+rotoidoscope")
            COMPREPLY=( $(compgen -W "> _center_x[%],_center_y[%],_tiles>0,_smoothness[%]>=0,_boundary_conditions={_0=dirichlet_|_1=neumann_|_2=periodic_|_3=mirror_}") ); return 0;;
        "round" | "-round" | "+round")
            COMPREPLY=( $(compgen -W "rounding_value>=0,_rounding_type (no_arg)") ); return 0;;
        "roundify" | "-roundify" | "+roundify")
            COMPREPLY=( $(compgen -W "> gamma>=0") ); return 0;;
        "rows" | "-rows" | "+rows")
            COMPREPLY=( $(compgen -W "> y0[%],_y1[%]") ); return 0;;
        "rprogress" | "-rprogress" | "+rprogress")
            COMPREPLY=( $(compgen -W "> 0<=value<=100_|_-1_|_\"command\",0<=value_min<=100,0<=value_max<=100") ); return 0;;
        "rr2d" | "-rr2d" | "+rr2d")
            COMPREPLY=( $(compgen -W "> width>0,height>0,_mode={_0=inside_|_1=outside_|_2=padded_},0=<_interpolation<=6") ); return 0;;
        "run" | "-run" | "+run")
            COMPREPLY=( $(compgen -W "> \"G\'MIC_pipeline\"") ); return 0;;
        "s3d" | "-s3d" | "+s3d")
            COMPREPLY=( $(compgen -W "> _full_split={_0_|_1_}") ); return 0;;
        "s" | "-s" | "+s")
            COMPREPLY=( $(compgen -W "{_x_|_y_|_z_|_c_}...{_x_|_y_|_z_|_c_},_split_mode keep_splitting_values={_+_|_-_},_{_x_|_y_|_z_|_c_}...{_x_|_y_|_z_|_c_},value1,_value2,... (no_arg)") ); return 0;;
        "sample" | "-sample" | "+sample")
            COMPREPLY=( $(compgen -W "_name1={_?_|_apples_|_balloons_|_barbara_|_boats_|_bottles_|_butterfly_|_cameraman_|_car_|_cat_|_cliff_|_chick_|_colorful_|_david_|_dog_|_duck_|_eagle_|_elephant_|_earth_|_flower_|_fruits_|_gmicky_|_gmicky_mahvin_|_gmicky_wilber_|_greece_|_gummy_|_house_|_inside_|_landscape_|_leaf_|_lena_|_leno_|_lion_|_mandrill_|_monalisa_|_monkey_|_parrots_|_pencils_|_peppers_|_portrait0_|_portrait1_|_portrait2_|_portrait3_|_portrait4_|_portrait5_|_portrait6_|_portrait7_|_portrait8_|_portrait9_|_roddy_|_rooster_|_rose_|_square_|_swan_|_teddy_|_tiger_|_tulips_|_wall_|_waterfall_|_zelda_},_name2,...,_nameN,_width={_>=0_|_0_(auto)_},_height_=_{_>=0_|_0_(auto)_} (no_arg)") ); return 0;;
        "scale_dcci2x" | "-scale_dcci2x" | "+scale_dcci2x")
            COMPREPLY=( $(compgen -W "> _edge_threshold>=0,_exponent>0,_extend_1px={_0=false_|_1=true_}") ); return 0;;
        "scanlines" | "-scanlines" | "+scanlines")
            COMPREPLY=( $(compgen -W "> _amplitude,_bandwidth,_shape={_0=block_|_1=triangle_|_2=sine_|_3=sine+_|_4=random_},_angle,_offset") ); return 0;;
        "screen" | "-screen" | "+screen")
            COMPREPLY=( $(compgen -W "> _x0[%],_y0[%],_x1[%],_y1[%]") ); return 0;;
        "seamcarve" | "-seamcarve" | "+seamcarve")
            COMPREPLY=( $(compgen -W "> _width[%]>=0,_height[%]>=0,_is_priority_channel={_0_|_1_},_is_antialiasing={_0_|_1_},_maximum_seams[%]>=0") ); return 0;;
        "segment_watershed" | "-segment_watershed" | "+segment_watershed")
            COMPREPLY=( $(compgen -W "> _threshold>=0") ); return 0;;
        "select" | "-select" | "+select")
            COMPREPLY=( $(compgen -W "> feature_type,_X[%]>=0,_Y[%]>=0,_Z[%]>=0,_exit_on_anykey={_0_|_1_},_is_deep_selection={_0_|_1_}") ); return 0;;
        "select_color" | "-select_color" | "+select_color")
            COMPREPLY=( $(compgen -W "> tolerance[%]>=0,col1,...,colN") ); return 0;;
        "serialize" | "-serialize" | "+serialize")
            COMPREPLY=( $(compgen -W "> _datatype,_is_compressed={_0_|_1_},_store_names={_0_|_1_}") ); return 0;;
        "set" | "-set" | "+set")
            COMPREPLY=( $(compgen -W "> value,_x[%],_y[%],_z[%],_c[%]") ); return 0;;
        "sh" | "-sh" | "+sh")
            COMPREPLY=( $(compgen -W "x0[%],x1[%],y[%],z[%],c[%] y0[%],y1[%],z[%],c[%] z0[%],z1[%],c[%] c0[%],c1[%] c0[%] (no_arg)") ); return 0;;
        "shade_stripes" | "-shade_stripes" | "+shade_stripes")
            COMPREPLY=( $(compgen -W "> _frequency>=0,_direction={_0=horizontal_|_1=vertical_},_darkness>=0,_lightness>=0") ); return 0;;
        "shadow_patch" | "-shadow_patch" | "+shadow_patch")
            COMPREPLY=( $(compgen -W "> _opacity>=0") ); return 0;;
        "shape2bump" | "-shape2bump" | "+shape2bump")
            COMPREPLY=( $(compgen -W "> _resolution>=0,0<=_weight_std_max_avg<=1,_dilation,_smoothness>=0") ); return 0;;
        "shape_circle" | "-shape_circle" | "+shape_circle")
            COMPREPLY=( $(compgen -W "> _size>=0") ); return 0;;
        "shape_cupid" | "-shape_cupid" | "+shape_cupid")
            COMPREPLY=( $(compgen -W "> _size>=0") ); return 0;;
        "shape_diamond" | "-shape_diamond" | "+shape_diamond")
            COMPREPLY=( $(compgen -W "> _size>=0") ); return 0;;
        "shape_dragon" | "-shape_dragon" | "+shape_dragon")
            COMPREPLY=( $(compgen -W "> _size>=0,_recursion_level>=0,_angle") ); return 0;;
        "shape_fern" | "-shape_fern" | "+shape_fern")
            COMPREPLY=( $(compgen -W "> _size>=0,_density[%]>=0,_angle,0<=_opacity<=1,_type={_0=Asplenium_adiantum-nigrum_|_1=Thelypteridaceae_}") ); return 0;;
        "shape_gear" | "-shape_gear" | "+shape_gear")
            COMPREPLY=( $(compgen -W "> _size>=0,_nb_teeth>0,0<=_height_teeth<=100,0<=_offset_teeth<=100,0<=_inner_radius<=100") ); return 0;;
        "shape_heart" | "-shape_heart" | "+shape_heart")
            COMPREPLY=( $(compgen -W "> _size>=0") ); return 0;;
        "shape_polygon" | "-shape_polygon" | "+shape_polygon")
            COMPREPLY=( $(compgen -W "> _size>=0,_nb_vertices>=3,_angle") ); return 0;;
        "shape_snowflake" | "-shape_snowflake" | "+shape_snowflake")
            COMPREPLY=( $(compgen -W "> size>=0,0<=_nb_recursions<=6") ); return 0;;
        "shape_star" | "-shape_star" | "+shape_star")
            COMPREPLY=( $(compgen -W "> _size>=0,_nb_branches>0,0<=_thickness<=1") ); return 0;;
        "shared" | "-shared" | "+shared")
            COMPREPLY=( $(compgen -W "x0[%],x1[%],y[%],z[%],c[%] y0[%],y1[%],z[%],c[%] z0[%],z1[%],c[%] c0[%],c1[%] c0[%] (no_arg)") ); return 0;;
        "sharpen" | "-sharpen" | "+sharpen")
            COMPREPLY=( $(compgen -W "amplitude>=0 amplitude>=0,edge>=0,_alpha[%],_sigma[%]") ); return 0;;
        "shift" | "-shift" | "+shift")
            COMPREPLY=( $(compgen -W "> vx[%],_vy[%],_vz[%],_vc[%],_boundary_conditions,_interpolation={_0=nearest_neighbor_|_1=linear_}") ); return 0;;
        "shift_tiles" | "-shift_tiles" | "+shift_tiles")
            COMPREPLY=( $(compgen -W "> M>0,_N>0,_amplitude") ); return 0;;
        "shrink_x" | "-shrink_x" | "+shrink_x")
            COMPREPLY=( $(compgen -W "> size_x>=0") ); return 0;;
        "shrink_xy" | "-shrink_xy" | "+shrink_xy")
            COMPREPLY=( $(compgen -W "> size>=0") ); return 0;;
        "shrink_xyz" | "-shrink_xyz" | "+shrink_xyz")
            COMPREPLY=( $(compgen -W "> size>=0") ); return 0;;
        "shrink_y" | "-shrink_y" | "+shrink_y")
            COMPREPLY=( $(compgen -W "> size_y>=0") ); return 0;;
        "shrink_z" | "-shrink_z" | "+shrink_z")
            COMPREPLY=( $(compgen -W "> size_z>=0") ); return 0;;
        "sierpinski3d" | "-sierpinski3d" | "+sierpinski3d")
            COMPREPLY=( $(compgen -W "> _recursion_level>=0,_width,_height") ); return 0;;
        "sierpinski" | "-sierpinski" | "+sierpinski")
            COMPREPLY=( $(compgen -W "> recursion_level>=0") ); return 0;;
        "skeleton3d" | "-skeleton3d" | "+skeleton3d")
            COMPREPLY=( $(compgen -W "> _metric,_frame_type={_0=squares_|_1=diamonds_|_2=circles_|_3=auto_},_skeleton_opacity,_frame_opacity,_is_frame_wireframe={_0_|_1_}") ); return 0;;
        "skeleton" | "-skeleton" | "+skeleton")
            COMPREPLY=( $(compgen -W "> _boundary_conditions={_0=dirichlet_|_1=neumann_}") ); return 0;;
        "sketchbw" | "-sketchbw" | "+sketchbw")
            COMPREPLY=( $(compgen -W "> _nb_angles>0,_start_angle,_angle_range>=0,_length>=0,_threshold>=0,_opacity,_bgfactor>=0,_density>0,_sharpness>=0,_anisotropy>=0,_smoothness>=0,_coherence>=0,_is_boost={_0_|_1_},_is_curved={_0_|_1_}") ); return 0;;
        "skip" | "-skip" | "+skip")
            COMPREPLY=( $(compgen -W "> item") ); return 0;;
        "sl3d" | "-sl3d" | "+sl3d")
            COMPREPLY=( $(compgen -W "> value>=0") ); return 0;;
        "slic" | "-slic" | "+slic")
            COMPREPLY=( $(compgen -W "> size>0,_regularity>=0,_nb_iterations>0") ); return 0;;
        "slices" | "-slices" | "+slices")
            COMPREPLY=( $(compgen -W "> z0[%],_z1[%]") ); return 0;;
        "smooth" | "-smooth" | "+smooth")
            COMPREPLY=( $(compgen -W "amplitude[%]>=0,_sharpness>=0,0<=_anisotropy<=1,_alpha[%],_sigma[%],_dl>0,_da>0,_precision>0,_interpolation,_fast_approx={_0_|_1_} nb_iterations>=0,_sharpness>=0,_anisotropy,_alpha,_sigma,_dt>0,0 [tensor_field],_amplitude>=0,_dl>0,_da>0,_precision>0,_interpolation,_fast_approx={_0_|_1_} [tensor_field],_nb_iters>=0,_dt>0,0") ); return 0;;
        "snapshot3d" | "-snapshot3d" | "+snapshot3d")
            COMPREPLY=( $(compgen -W "_size>0,_zoom>=0,_backgroundR,_backgroundG,_backgroundB,_backgroundA [background_image],zoom>=0") ); return 0;;
        "solidify" | "-solidify" | "+solidify")
            COMPREPLY=( $(compgen -W "> _smoothness[%]>=0,_diffusion_type={_0=isotropic_|_1=Delaunay-guided_|_2=edge-oriented_},_diffusion_iter>=0") ); return 0;;
        "solve" | "-solve" | "+solve")
            COMPREPLY=( $(compgen -W "> [image]") ); return 0;;
        "solve_poisson" | "-solve_poisson" | "+solve_poisson")
            COMPREPLY=( $(compgen -W "> \"laplacian_command\",_nb_iterations>=0,_time_step>0,_nb_scales>=0") ); return 0;;
        "sort" | "-sort" | "+sort")
            COMPREPLY=( $(compgen -W "> _ordering={_+_|_-_},_axis={_x_|_y_|_z_|_c_}") ); return 0;;
        "sort_list" | "-sort_list" | "+sort_list")
            COMPREPLY=( $(compgen -W "> _ordering={_+_|_-_},_criterion") ); return 0;;
        "sp" | "-sp" | "+sp")
            COMPREPLY=( $(compgen -W "_name1={_?_|_apples_|_balloons_|_barbara_|_boats_|_bottles_|_butterfly_|_cameraman_|_car_|_cat_|_cliff_|_chick_|_colorful_|_david_|_dog_|_duck_|_eagle_|_elephant_|_earth_|_flower_|_fruits_|_gmicky_|_gmicky_mahvin_|_gmicky_wilber_|_greece_|_gummy_|_house_|_inside_|_landscape_|_leaf_|_lena_|_leno_|_lion_|_mandrill_|_monalisa_|_monkey_|_parrots_|_pencils_|_peppers_|_portrait0_|_portrait1_|_portrait2_|_portrait3_|_portrait4_|_portrait5_|_portrait6_|_portrait7_|_portrait8_|_portrait9_|_roddy_|_rooster_|_rose_|_square_|_swan_|_teddy_|_tiger_|_tulips_|_wall_|_waterfall_|_zelda_},_name2,...,_nameN,_width={_>=0_|_0_(auto)_},_height_=_{_>=0_|_0_(auto)_} (no_arg)") ); return 0;;
        "specl3d" | "-specl3d" | "+specl3d")
            COMPREPLY=( $(compgen -W "> value>=0") ); return 0;;
        "specs3d" | "-specs3d" | "+specs3d")
            COMPREPLY=( $(compgen -W "> value>=0") ); return 0;;
        "sphere3d" | "-sphere3d" | "+sphere3d")
            COMPREPLY=( $(compgen -W "> radius,_nb_recursions>=0") ); return 0;;
        "spherical3d" | "-spherical3d" | "+spherical3d")
            COMPREPLY=( $(compgen -W "> _nb_azimuth>=3,_nb_zenith>=3,_radius_function(phi,theta)") ); return 0;;
        "spherize" | "-spherize" | "+spherize")
            COMPREPLY=( $(compgen -W "> _radius[%]>=0,_strength,_smoothness[%]>=0,_center_x[%],_center_y[%],_ratio_x/y>0,_angle,_interpolation") ); return 0;;
        "spiralbw" | "-spiralbw" | "+spiralbw")
            COMPREPLY=( $(compgen -W "> width>0,_height>0,_is_2dcoords={_0_|_1_}") ); return 0;;
        "spline3d" | "-spline3d" | "+spline3d")
            COMPREPLY=( $(compgen -W "> x0[%],y0[%],z0[%],u0[%],v0[%],w0[%],x1[%],y1[%],z1[%],u1[%],v1[%],w1[%],_nb_vertices>=2") ); return 0;;
        "spline" | "-spline" | "+spline")
            COMPREPLY=( $(compgen -W "> x0[%],y0[%],u0[%],v0[%],x1[%],y1[%],u1[%],v1[%],_opacity,_color1,...") ); return 0;;
        "split3d" | "-split3d" | "+split3d")
            COMPREPLY=( $(compgen -W "> _full_split={_0_|_1_}") ); return 0;;
        "split" | "-split" | "+split")
            COMPREPLY=( $(compgen -W "{_x_|_y_|_z_|_c_}...{_x_|_y_|_z_|_c_},_split_mode keep_splitting_values={_+_|_-_},_{_x_|_y_|_z_|_c_}...{_x_|_y_|_z_|_c_},value1,_value2,... (no_arg)") ); return 0;;
        "split_colors" | "-split_colors" | "+split_colors")
            COMPREPLY=( $(compgen -W "> _tolerance>=0,_max_nb_outputs>0,_min_area>0") ); return 0;;
        "split_details" | "-split_details" | "+split_details")
            COMPREPLY=( $(compgen -W "> _nb_scales>0,_base_scale[%]>=0,_detail_scale[%]>=0") ); return 0;;
        "split_freq" | "-split_freq" | "+split_freq")
            COMPREPLY=( $(compgen -W "> smoothness>0[%]") ); return 0;;
        "split_tiles" | "-split_tiles" | "+split_tiles")
            COMPREPLY=( $(compgen -W "> M!=0,_N!=0,_is_homogeneous={_0_|_1_}") ); return 0;;
        "sponge" | "-sponge" | "+sponge")
            COMPREPLY=( $(compgen -W "> _size>0") ); return 0;;
        "spread" | "-spread" | "+spread")
            COMPREPLY=( $(compgen -W "> _dx>=0,_dy>=0,_dz>=0") ); return 0;;
        "sprites3d" | "-sprites3d" | "+sprites3d")
            COMPREPLY=( $(compgen -W "> [sprite],_sprite_has_alpha_channel={_0_|_1_}") ); return 0;;
        "srand" | "-srand" | "+srand")
            COMPREPLY=( $(compgen -W "value (no_arg)") ); return 0;;
        "srgb2lab8" | "-srgb2lab8" | "+srgb2lab8")
            COMPREPLY=( $(compgen -W "illuminant={_0=D50_|_1=D65_|_2=E_} (no_arg)") ); return 0;;
        "srgb2lab" | "-srgb2lab" | "+srgb2lab")
            COMPREPLY=( $(compgen -W "illuminant={_0=D50_|_1=D65_|_2=E_} (no_arg)") ); return 0;;
        "ss3d" | "-ss3d" | "+ss3d")
            COMPREPLY=( $(compgen -W "> value>=0") ); return 0;;
        "ssd_patch" | "-ssd_patch" | "+ssd_patch")
            COMPREPLY=( $(compgen -W "> [patch],_use_fourier={_0_|_1_},_boundary_conditions") ); return 0;;
        "ssim" | "-ssim" | "+ssim")
            COMPREPLY=( $(compgen -W "> [reference],_patch_size>0,_max_value>0") ); return 0;;
        "ssim_matrix" | "-ssim_matrix" | "+ssim_matrix")
            COMPREPLY=( $(compgen -W "> _patch_size>0,_max_value>0") ); return 0;;
        "stained_glass" | "-stained_glass" | "+stained_glass")
            COMPREPLY=( $(compgen -W "> _edges[%]>=0,_shading>=0,_is_thin_separators={_0_|_1_}") ); return 0;;
        "star3d" | "-star3d" | "+star3d")
            COMPREPLY=( $(compgen -W "> _nb_branches>0,0<=_thickness<=1") ); return 0;;
        "stars" | "-stars" | "+stars")
            COMPREPLY=( $(compgen -W "> _density[%]>=0,_depth>=0,_size>0,_nb_branches>=1,0<=_thickness<=1,_smoothness[%]>=0,_R,_G,_B,_opacity") ); return 0;;
        "status" | "-status" | "+status")
            COMPREPLY=( $(compgen -W "> status_string") ); return 0;;
        "stencil" | "-stencil" | "+stencil")
            COMPREPLY=( $(compgen -W "> _radius[%]>=0,_smoothness>=0,_iterations>=0") ); return 0;;
        "stencilbw" | "-stencilbw" | "+stencilbw")
            COMPREPLY=( $(compgen -W "> _edges>=0,_smoothness>=0") ); return 0;;
        "store" | "-store" | "+store")
            COMPREPLY=( $(compgen -W "> _is_compressed={_0_|_1_},variable_name1,_variable_name2,...") ); return 0;;
        "str2hex" | "-str2hex" | "+str2hex")
            COMPREPLY=( $(compgen -W "> \"string\"") ); return 0;;
        "str" | "-str" | "+str")
            COMPREPLY=( $(compgen -W "> string") ); return 0;;
        "strcapitalize" | "-strcapitalize" | "+strcapitalize")
            COMPREPLY=( $(compgen -W "> string") ); return 0;;
        "strcasevar" | "-strcasevar" | "+strcasevar")
            COMPREPLY=( $(compgen -W "> \"string\"") ); return 0;;
        "strcontains" | "-strcontains" | "+strcontains")
            COMPREPLY=( $(compgen -W "> string1,string2") ); return 0;;
        "streamline3d" | "-streamline3d" | "+streamline3d")
            COMPREPLY=( $(compgen -W "x[%],y[%],z[%],_L>=0,_dl>0,_interpolation,_is_backward={_0_|_1_},_is_oriented={_0_|_1_} \'formula\',x,y,z,_L>=0,_dl>0,_interpolation,_is_backward={_0_|_1_},_is_oriented={_0_|_1_}") ); return 0;;
        "stripes_y" | "-stripes_y" | "+stripes_y")
            COMPREPLY=( $(compgen -W "> _frequency>=0") ); return 0;;
        "strlen" | "-strlen" | "+strlen")
            COMPREPLY=( $(compgen -W "> string1") ); return 0;;
        "strlowercase" | "-strlowercase" | "+strlowercase")
            COMPREPLY=( $(compgen -W "> string") ); return 0;;
        "strreplace" | "-strreplace" | "+strreplace")
            COMPREPLY=( $(compgen -W "> string,search,replace") ); return 0;;
        "structuretensors" | "-structuretensors" | "+structuretensors")
            COMPREPLY=( $(compgen -W "> _scheme={_0=centered_|_1=forward/backward_}") ); return 0;;
        "struppercase" | "-struppercase" | "+struppercase")
            COMPREPLY=( $(compgen -W "> string") ); return 0;;
        "strvar" | "-strvar" | "+strvar")
            COMPREPLY=( $(compgen -W "> \"string\"") ); return 0;;
        "strver" | "-strver" | "+strver")
            COMPREPLY=( $(compgen -W "> _version,_prerelease") ); return 0;;
        "stylize" | "-stylize" | "+stylize")
            COMPREPLY=( $(compgen -W "> [style_image],_fidelity_finest,_fidelity_coarsest,_fidelity_smoothness_finest>=0,_fidelity_smoothnes_coarsest>=0,0<=_fidelity_chroma<=1,_init_type,_init_resolution>=0,init_max_gradient>=0,_patchsize_analysis>0,_patchsize_synthesis>0,_patchsize_synthesis_final>0,_nb_matches_finest>=0,_nb_matches_coarsest>=0,_penalize_repetitions>=0,_matching_precision>=0,_scale_factor>1,_skip_finest_scales>=0,_\"image_matching_command\"") ); return 0;;
        "sub3d" | "-sub3d" | "+sub3d")
            COMPREPLY=( $(compgen -W "> tx,_ty,_tz") ); return 0;;
        "sub" | "-sub" | "+sub")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "sub_alpha" | "-sub_alpha" | "+sub_alpha")
            COMPREPLY=( $(compgen -W "> [base_image],_opacity_gain>=1") ); return 0;;
        "superformula3d" | "-superformula3d" | "+superformula3d")
            COMPREPLY=( $(compgen -W "> resolution>1,m>=1,n1,n2,n3") ); return 0;;
        "surfels3d" | "-surfels3d" | "+surfels3d")
            COMPREPLY=( $(compgen -W "> 0<=_left_right_attenuation<=1,0<=_top_bottom_attenuation<=1,0<=_closer_further_attenuation<=1") ); return 0;;
        "symmetrize" | "-symmetrize" | "+symmetrize")
            COMPREPLY=( $(compgen -W "> _x[%],_y[%],_angle,_boundary_conditions={_0=dirichlet_|_1=neumann_|_2=periodic_|_3=mirror_},_is_antisymmetry={_0_|_1_},_swap_sides={_0_|_1_}") ); return 0;;
        "syntexturize" | "-syntexturize" | "+syntexturize")
            COMPREPLY=( $(compgen -W "> _width[%]>0,_height[%]>0") ); return 0;;
        "syntexturize_matchpatch" | "-syntexturize_matchpatch" | "+syntexturize_matchpatch")
            COMPREPLY=( $(compgen -W "> _width[%]>0,_height[%]>0,_nb_scales>=0,_patch_size>0,_blending_size>=0,_precision>=0") ); return 0;;
        "t3d" | "-t3d" | "+t3d")
            COMPREPLY=( $(compgen -W "> [ind_texture],_[ind_coords]") ); return 0;;
        "t" | "-t" | "+t")
            COMPREPLY=( $(compgen -W "> text,_x[%|~],_y[%|~],_font_height[%]>=0,_opacity,_color1,...") ); return 0;;
        "taquin" | "-taquin" | "+taquin")
            COMPREPLY=( $(compgen -W "> M>0,_N>0,_remove_tile={_0=none_|_1=first_|_2=last_|_3=random_},_relief,_border_thickness[%],_border_outline[%],_outline_color") ); return 0;;
        "tensors3d" | "-tensors3d" | "+tensors3d")
            COMPREPLY=( $(compgen -W "> _radius_factor>=0,_shape={_0=box_|_>=N=ellipsoid_},_radius_min>=0") ); return 0;;
        "testimage2d" | "-testimage2d" | "+testimage2d")
            COMPREPLY=( $(compgen -W "> _width>0,_height>0,_spectrum>0") ); return 0;;
        "tetraedron_shade" | "-tetraedron_shade" | "+tetraedron_shade")
            COMPREPLY=( $(compgen -W "> x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3,R0,G0,B0,...,R1,G1,B1,...,R2,G2,B2,...,R3,G3,B3,...") ); return 0;;
        "tetris" | "-tetris" | "+tetris")
            COMPREPLY=( $(compgen -W "> _scale>0") ); return 0;;
        "text3d" | "-text3d" | "+text3d")
            COMPREPLY=( $(compgen -W "> text,_font_height>0,_depth>0,_smoothness") ); return 0;;
        "text" | "-text" | "+text")
            COMPREPLY=( $(compgen -W "> text,_x[%|~],_y[%|~],_font_height[%]>=0,_opacity,_color1,...") ); return 0;;
        "text_outline" | "-text_outline" | "+text_outline")
            COMPREPLY=( $(compgen -W "> text,_x[%|~],_y[%|~],_font_height[%]>0,_outline>=0,_opacity,_color1,...") ); return 0;;
        "text_pointcloud3d" | "-text_pointcloud3d" | "+text_pointcloud3d")
            COMPREPLY=( $(compgen -W "> _\"text1\",_\"text2\",_smoothness") ); return 0;;
        "texturize3d" | "-texturize3d" | "+texturize3d")
            COMPREPLY=( $(compgen -W "> [ind_texture],_[ind_coords]") ); return 0;;
        "texturize_canvas" | "-texturize_canvas" | "+texturize_canvas")
            COMPREPLY=( $(compgen -W "> _amplitude>=0,_fibrousness>=0,_emboss_level>=0") ); return 0;;
        "thickline" | "-thickline" | "+thickline")
            COMPREPLY=( $(compgen -W "> x0[%],y0[%],x1[%],y1[%],_thickness,_opacity,_color1") ); return 0;;
        "thinning" | "-thinning" | "+thinning")
            COMPREPLY=( $(compgen -W "> _boundary_conditions={_0=dirichlet_|_1=neumann_}") ); return 0;;
        "threshold" | "-threshold" | "+threshold")
            COMPREPLY=( $(compgen -W "> value[%],_is_soft={_0_|_1_}") ); return 0;;
        "tixy" | "-tixy" | "+tixy")
            COMPREPLY=( $(compgen -W "> \"expression\"") ); return 0;;
        "to" | "-to" | "+to")
            COMPREPLY=( $(compgen -W "> text,_x[%|~],_y[%|~],_font_height[%]>0,_outline>=0,_opacity,_color1,...") ); return 0;;
        "to_clutname" | "-to_clutname" | "+to_clutname")
            COMPREPLY=( $(compgen -W "> \"string\"") ); return 0;;
        "to_colormode" | "-to_colormode" | "+to_colormode")
            COMPREPLY=( $(compgen -W "> mode={_0=adaptive_|_1=G_|_2=GA_|_3=RGB_|_4=RGBA_}") ); return 0;;
        "to_pseudogray" | "-to_pseudogray" | "+to_pseudogray")
            COMPREPLY=( $(compgen -W "> _max_step>=0,_is_perceptual_constraint={_0_|_1_},_bits_depth>0") ); return 0;;
        "tones" | "-tones" | "+tones")
            COMPREPLY=( $(compgen -W "> N>0") ); return 0;;
        "topographic_map" | "-topographic_map" | "+topographic_map")
            COMPREPLY=( $(compgen -W "> _nb_levels>0,_smoothness") ); return 0;;
        "torus3d" | "-torus3d" | "+torus3d")
            COMPREPLY=( $(compgen -W "> _radius1,_radius2,_nb_subdivisions1>2,_nb_subdivisions2>2") ); return 0;;
        "transfer_histogram" | "-transfer_histogram" | "+transfer_histogram")
            COMPREPLY=( $(compgen -W "> [reference_image],_nb_levels>0,_color_channels") ); return 0;;
        "transfer_pca" | "-transfer_pca" | "+transfer_pca")
            COMPREPLY=( $(compgen -W "> [reference_image],_color_channels") ); return 0;;
        "transfer_rgb" | "-transfer_rgb" | "+transfer_rgb")
            COMPREPLY=( $(compgen -W "> [target],_gamma>=0,_regularization>=0,_luminosity_constraints>=0,_rgb_resolution>=0,_is_constraints={_0_|_1_}") ); return 0;;
        "transform_polar" | "-transform_polar" | "+transform_polar")
            COMPREPLY=( $(compgen -W "> \"expr_radius\",_\"expr_angle\",_center_x[%],_center_y[%],_boundary_conditions={_0=dirichlet_|_1=neumann_|_2=periodic_|_3=mirror_}") ); return 0;;
        "transition3d" | "-transition3d" | "+transition3d")
            COMPREPLY=( $(compgen -W "> _nb_frames>=2,_nb_xtiles>0,_nb_ytiles>0,_axis_x,_axis_y,_axis_z,_is_antialias={_0_|_1_}") ); return 0;;
        "transition" | "-transition" | "+transition")
            COMPREPLY=( $(compgen -W "> [transition_shape],nb_added_frames>=0,100>=shading>=0,_single_frame_only={_-1=disabled_|_>=0_}") ); return 0;;
        "triangle3d" | "-triangle3d" | "+triangle3d")
            COMPREPLY=( $(compgen -W "> x0,y0,z0,x1,y1,z1,x2,y2,z2") ); return 0;;
        "triangle_shade" | "-triangle_shade" | "+triangle_shade")
            COMPREPLY=( $(compgen -W "> x0,y0,x1,y1,x2,y2,R0,G0,B0,...,R1,G1,B1,...,R2,G2,B2,...") ); return 0;;
        "trisolve" | "-trisolve" | "+trisolve")
            COMPREPLY=( $(compgen -W "> [image]") ); return 0;;
        "truchet" | "-truchet" | "+truchet")
            COMPREPLY=( $(compgen -W "> _scale>0,_radius>=0,_pattern_type={_0=straight_|_1=curved_}") ); return 0;;
        "tsp" | "-tsp" | "+tsp")
            COMPREPLY=( $(compgen -W "> _precision>=0") ); return 0;;
        "tunnel" | "-tunnel" | "+tunnel")
            COMPREPLY=( $(compgen -W "> _level>=0,_factor>0,_centering_x,_centering_y,_opacity,_angle") ); return 0;;
        "turbulence" | "-turbulence" | "+turbulence")
            COMPREPLY=( $(compgen -W "> _radius>0,_octaves={1,2,3...,12},_alpha>0,_difference={-10,10},_mode={0,1,2,3}") ); return 0;;
        "tv_flow" | "-tv_flow" | "+tv_flow")
            COMPREPLY=( $(compgen -W "> _nb_iter>=0,_dt,_keep_sequence={_0_|_1_}") ); return 0;;
        "twirl" | "-twirl" | "+twirl")
            COMPREPLY=( $(compgen -W "> _amplitude,_center_x[%],_center_y[%],_boundary_conditions={_0=dirichlet_|_1=neumann_|_2=periodic_|_3=mirror_}") ); return 0;;
        "u" | "-u" | "+u")
            COMPREPLY=( $(compgen -W "> status_string") ); return 0;;
        "uint82base64" | "-uint82base64" | "+uint82base64")
            COMPREPLY=( $(compgen -W "> _encoding={_0=base64_|_1=base64url_}") ); return 0;;
        "um" | "-um" | "+um")
            COMPREPLY=( $(compgen -W "command_name[,_command_name2,...] *") ); return 0;;
        "uncommand" | "-uncommand" | "+uncommand")
            COMPREPLY=( $(compgen -W "command_name[,_command_name2,...] *") ); return 0;;
        "undistort" | "-undistort" | "+undistort")
            COMPREPLY=( $(compgen -W "> -1<=_amplitude<=1,_aspect_ratio,_zoom,_center_x[%],_center_y[%],_boundary_conditions") ); return 0;;
        "uniform_distribution" | "-uniform_distribution" | "+uniform_distribution")
            COMPREPLY=( $(compgen -W "> nb_levels>=1,spectrum>=1") ); return 0;;
        "unroll" | "-unroll" | "+unroll")
            COMPREPLY=( $(compgen -W "> _axis={_x_|_y_|_z_|_c_}") ); return 0;;
        "unserialize" | "-unserialize" | "+unserialize")
            COMPREPLY=( $(compgen -W "> ") ); return 0;;
        "unsharp" | "-unsharp" | "+unsharp")
            COMPREPLY=( $(compgen -W "> radius[%]>=0,_amount>=0,_threshold[%]>=0") ); return 0;;
        "unsharp_octave" | "-unsharp_octave" | "+unsharp_octave")
            COMPREPLY=( $(compgen -W "> _nb_scales>0,_radius[%]>=0,_amount>=0,threshold[%]>=0") ); return 0;;
        "upscale_smart" | "-upscale_smart" | "+upscale_smart")
            COMPREPLY=( $(compgen -W "> width[%],_height[%],_depth,_smoothness>=0,_anisotropy=[0,1],sharpening>=0") ); return 0;;
        "v" | "-v" | "+v")
            COMPREPLY=( $(compgen -W "level {_+_|_-_}") ); return 0;;
        "vanvliet" | "-vanvliet" | "+vanvliet")
            COMPREPLY=( $(compgen -W "> std_deviation>=0[%],order={_0_|_1_|_2_|_3_},axis={_x_|_y_|_z_|_c_},_boundary_conditions") ); return 0;;
        "variance_patch" | "-variance_patch" | "+variance_patch")
            COMPREPLY=( $(compgen -W "> _patch_size>=1") ); return 0;;
        "verbose" | "-verbose" | "+verbose")
            COMPREPLY=( $(compgen -W "level {_+_|_-_}") ); return 0;;
        "video2files" | "-video2files" | "+video2files")
            COMPREPLY=( $(compgen -W "> input_filename,_output_filename,_first_frame>=0,_last_frame={_>=0_|_-1=last_},_frame_step>=1") ); return 0;;
        "vignette" | "-vignette" | "+vignette")
            COMPREPLY=( $(compgen -W "> _strength>=0,0<=_radius_min<=100,0<=_radius_max<=100") ); return 0;;
        "w" | "-w" | "+w")
            COMPREPLY=( $(compgen -W "> _width[%]>=-1,_height[%]>=-1,_normalization,_fullscreen,_pos_x[%],_pos_y[%],_title") ); return 0;;
        "wait" | "-wait" | "+wait")
            COMPREPLY=( $(compgen -W "delay (no_arg)") ); return 0;;
        "warhol" | "-warhol" | "+warhol")
            COMPREPLY=( $(compgen -W "> _M>0,_N>0,_smoothness>=0,_color>=0") ); return 0;;
        "warn" | "-warn" | "+warn")
            COMPREPLY=( $(compgen -W "> _force_visible={_0_|_1_},_message") ); return 0;;
        "warp" | "-warp" | "+warp")
            COMPREPLY=( $(compgen -W "> [warping_field],_mode,_interpolation,_boundary_conditions,_nb_frames>0") ); return 0;;
        "warp_patch" | "-warp_patch" | "+warp_patch")
            COMPREPLY=( $(compgen -W "> [warping_field],patch_width>=1,_patch_height>=1,_patch_depth>=1,_std_factor>0,_boundary_conditions.") ); return 0;;
        "warp_perspective" | "-warp_perspective" | "+warp_perspective")
            COMPREPLY=( $(compgen -W "> _x-angle,_y-angle,_zoom>0,_x-center,_y-center,_boundary_conditions={_0=dirichlet_|_1=neumann_|_2=periodic_|_3=mirror_}") ); return 0;;
        "warp_rbf" | "-warp_rbf" | "+warp_rbf")
            COMPREPLY=( $(compgen -W "> xs0[%],ys0[%],xt0[%],yt0[%],...,xsN[%],ysN[%],xtN[%],ytN[%]") ); return 0;;
        "water" | "-water" | "+water")
            COMPREPLY=( $(compgen -W "> _amplitude,_smoothness>=0,_angle") ); return 0;;
        "watermark_fourier" | "-watermark_fourier" | "+watermark_fourier")
            COMPREPLY=( $(compgen -W "> text,_size>0") ); return 0;;
        "watermark_visible" | "-watermark_visible" | "+watermark_visible")
            COMPREPLY=( $(compgen -W "> _text,0<_opacity<1,_size>0,_angle,_mode={_0=remove_|_1=add_},_smoothness>=0") ); return 0;;
        "watershed" | "-watershed" | "+watershed")
            COMPREPLY=( $(compgen -W "> [priority_image],_is_high_connectivity={_0_|_1_}") ); return 0;;
        "wave" | "-wave" | "+wave")
            COMPREPLY=( $(compgen -W "> _amplitude>=0,_frequency>=0,_center_x,_center_y") ); return 0;;
        "weave" | "-weave" | "+weave")
            COMPREPLY=( $(compgen -W "> _density>=0,0<=_thickness<=100,0<=_shadow<=100,_shading>=0,_fibers_amplitude>=0,_fibers_smoothness>=0,_angle,-1<=_x_curvature<=1,-1<=_y_curvature<=1") ); return 0;;
        "weird3d" | "-weird3d" | "+weird3d")
            COMPREPLY=( $(compgen -W "> _resolution>0") ); return 0;;
        "while" | "-while" | "+while")
            COMPREPLY=( $(compgen -W "> condition") ); return 0;;
        "whirls" | "-whirls" | "+whirls")
            COMPREPLY=( $(compgen -W "> _texture>=0,_smoothness>=0,_darkness>=0,_lightness>=0") ); return 0;;
        "wind" | "-wind" | "+wind")
            COMPREPLY=( $(compgen -W "> _amplitude>=0,_angle,0<=_attenuation<=1,_threshold") ); return 0;;
        "window" | "-window" | "+window")
            COMPREPLY=( $(compgen -W "> _width[%]>=-1,_height[%]>=-1,_normalization,_fullscreen,_pos_x[%],_pos_y[%],_title") ); return 0;;
        "x" | "-x" | "+x")
            COMPREPLY=( $(compgen -W "> _is_verbose={_0_|_1_},\"command\"") ); return 0;;
        "x_color_curves" | "-x_color_curves" | "+x_color_curves")
            COMPREPLY=( $(compgen -W "> _colorspace={_rgb_|_cmy_|_cmyk_|_hsi_|_hsl_|_hsv_|_lab_|_lch_|_ycbcr_|_last_}") ); return 0;;
        "x_colorize" | "-x_colorize" | "+x_colorize")
            COMPREPLY=( $(compgen -W "> _is_lineart={_0_|_1_},_max_resolution={_0_|_>=128_},_multichannels_output={_0_|_1_},_[palette1],_[palette2],_[grabber1]") ); return 0;;
        "x_grab_color" | "-x_grab_color" | "+x_grab_color")
            COMPREPLY=( $(compgen -W "> _variable_name") ); return 0;;
        "x_jawbreaker" | "-x_jawbreaker" | "+x_jawbreaker")
            COMPREPLY=( $(compgen -W "> 0<_width<20,0<_height<20,0<_balls<=8") ); return 0;;
        "x_mandelbrot" | "-x_mandelbrot" | "+x_mandelbrot")
            COMPREPLY=( $(compgen -W "> _julia={_0_|_1_},_c0r,_c0i") ); return 0;;
        "x_mask_color" | "-x_mask_color" | "+x_mask_color")
            COMPREPLY=( $(compgen -W "> _colorspace={_all_|_rgb_|_lrgb_|_ycbcr_|_lab_|_lch_|_hsv_|_hsi_|_hsl_|_cmy_|_cmyk_|_yiq_},_spatial_tolerance>=0,_color_tolerance>=0") ); return 0;;
        "x_minesweeper" | "-x_minesweeper" | "+x_minesweeper")
            COMPREPLY=( $(compgen -W "> 8<=_width=<20,8<=_height<=20") ); return 0;;
        "x_morph" | "-x_morph" | "+x_morph")
            COMPREPLY=( $(compgen -W "> _nb_frames>=2,_preview_fidelity={_0=coarsest_|_1=coarse_|_2=normal_|_3=fine_|_4=finest_}") ); return 0;;
        "x_quantize_rgb" | "-x_quantize_rgb" | "+x_quantize_rgb")
            COMPREPLY=( $(compgen -W "> _nbcolors>=2") ); return 0;;
        "x_segment" | "-x_segment" | "+x_segment")
            COMPREPLY=( $(compgen -W "> _max_resolution={_0_|_>=128_}") ); return 0;;
        "x_select_color" | "-x_select_color" | "+x_select_color")
            COMPREPLY=( $(compgen -W "> _variable_name") ); return 0;;
        "x_select_function1d" | "-x_select_function1d" | "+x_select_function1d")
            COMPREPLY=( $(compgen -W "> _variable_name,_background_curve_R,_background_curve_G,_background_curve_B") ); return 0;;
        "x_select_palette" | "-x_select_palette" | "+x_select_palette")
            COMPREPLY=( $(compgen -W "> _variable_name,_number_of_columns={_0=auto_|_>0_}") ); return 0;;
        "x_warp" | "-x_warp" | "+x_warp")
            COMPREPLY=( $(compgen -W "> _nb_keypoints_xgrid>=2,_nb_keypoints_ygrid>=2,_nb_keypoints_contours>=0,_preview_fidelity={_0=coarsest_|_1=coarse_|_2=normal_|_3=fine_|_4=finest_},_[background_image],0<=_background_opacity<=1") ); return 0;;
        "x_whirl" | "-x_whirl" | "+x_whirl")
            COMPREPLY=( $(compgen -W "> _opacity>=0") ); return 0;;
        "xo" | "-xo" | "+xo")
            COMPREPLY=( $(compgen -W "> _mode,\"command\"") ); return 0;;
        "xor" | "-xor" | "+xor")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
        "xyz2lab" | "-xyz2lab" | "+xyz2lab")
            COMPREPLY=( $(compgen -W "illuminant={_0=D50_|_1=D65_|_2=E_} (no_arg)") ); return 0;;
        "xyz2rgb" | "-xyz2rgb" | "+xyz2rgb")
            COMPREPLY=( $(compgen -W "illuminant={_0=D50_|_1=D65_|_2=E_} (no_arg)") ); return 0;;
        "xyz82rgb" | "-xyz82rgb" | "+xyz82rgb")
            COMPREPLY=( $(compgen -W "illuminant={_0=D50_|_1=D65_|_2=E_} (no_arg)") ); return 0;;
        "y" | "-y" | "+y")
            COMPREPLY=( $(compgen -W "> _axis={_x_|_y_|_z_|_c_}") ); return 0;;
        "z" | "-z" | "+z")
            COMPREPLY=( $(compgen -W "x0[%],x1[%],_boundary_conditions x0[%],y0[%],x1[%],y1[%],_boundary_conditions x0[%],y0[%],z0[%],x1[%],y1[%],z1[%],_boundary_conditions x0[%],y0[%],z0[%],c0[%],x1[%],y1[%],z1[%],c1[%],_boundary_conditions") ); return 0;;
        "zoom" | "-zoom" | "+zoom")
            COMPREPLY=( $(compgen -W "> _factor,_cx,_cy,_cz,_boundary_conditions={_0=dirichlet_|_1=neumann_|_2=periodic_|_3=mirror_}") ); return 0;;
        "|" | "-|" | "+|")
            COMPREPLY=( $(compgen -W "value[%] [image] \'formula\' (no_arg)") ); return 0;;
    esac

    COMPREPLY=( $(compgen -W "$opts" -- "$cur") )
    if type -t _filedir >/dev/null; then
        _filedir
    else
        comptopt -o filenames 2>/dev/null
        COMPREPLY=( $(compgen -f -- ${cur}) )
    fi
}
complete -F _gmic -o filenames gmic