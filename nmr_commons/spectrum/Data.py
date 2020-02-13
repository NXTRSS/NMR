import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm


class Data(object):
    #
    # Constructor
    def __init__(self):
        self.data = None

    @staticmethod
    def visualize_data(data, output_path=None, points_x=None, points_y=None, points_x_2=None, points_y_2=None,
                       points_x_3=None, points_y_3=None, first_contour=None, number_of_contours=30, title=None, axes_obj=None):

        if axes_obj is None:
            plt.figure()

        if title is not None:
            plt.title(title, fontsize=10)

        if first_contour is None:
            first_contour = np.median(data)

        first_contour = max(0.0001, first_contour)

        levels_m = list(map(lambda x: -first_contour * (1.4 ** x), range(number_of_contours, 1, -1)))
        levels = list(map(lambda x: first_contour * (1.4 ** x), range(1, number_of_contours, 1)))

        if axes_obj is None:
            plt.contour(data, levels_m + levels, linewidths=1, cmap=cm.gray)
        else:
            axes_obj.contour(data, levels_m + levels, linewidths=1, cmap=cm.gray)

        plt.hold(True)

        if points_x is not None:

            if axes_obj is None:
                plt.scatter(points_x, points_y, marker='o', c='r', s=10, zorder=10, lw=0)
            else:
                axes_obj.scatter(points_x, points_y, marker='o', c='r', s=10, zorder=10, lw=0)

        if points_x_2 is not None:

            if axes_obj is None:
                plt.scatter(points_x_2, points_y_2, marker='o', c='b', s=10, zorder=10, lw=0)
            else:
                axes_obj.scatter(points_x_2, points_y_2, marker='o', c='b', s=10, zorder=10, lw=0)

        if points_x_3 is not None:

            if axes_obj is None:
                plt.scatter(points_x_3, points_y_3, marker='o', c='g', s=10, zorder=10, lw=0)
            else:
                axes_obj.scatter(points_x_3, points_y_3, marker='o', c='g', s=10, zorder=10, lw=0)

        if output_path is not None:
            manager = plt.get_current_fig_manager()
            manager.window.showMaximized()
            fig = plt.gcf()
            fig.set_size_inches(24, 14)
            fig.savefig(output_path, dpi=200)
