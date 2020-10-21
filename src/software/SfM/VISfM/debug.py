import os
import LayoutBuilder

def vis_sfm_json(sfm_json_bin_path):
    if not os.path.exists(sfm_json_bin_path):
        print('%s not exist' % sfm_json_bin_path)
        return

    dir, ext = os.path.splitext(sfm_json_bin_path)
    pc_dir = os.path.join(dir, 'point_cloud')

    ws = LayoutBuilder.WorkSpace(dir)

    if not os.path.exists(pc_dir):
        ws.openMVG_main_openMVG2Colmap(sfm_json_bin_path, dir)
        ws.colmap_model_converter(dir, dir, 'BIN')
        ws.colmap_2_point_cloud(dir, pc_dir)

    ws.vis_point_cloud(dir)

vis_sfm_json(
    'D:/xin/test/sfm_data.bin'
)

vis_sfm_json(
    'D:/xin/result/Allresutl.bin'
)