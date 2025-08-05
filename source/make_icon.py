#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon
import matplotlib.patches as patches

def create_hexagon(center_x=0, center_y=0, radius=1):
    """Create hexagon vertices"""
    angles = np.linspace(0, 2*np.pi, 7)
    x = center_x + radius * np.cos(angles)
    y = center_y + radius * np.sin(angles)
    return list(zip(x, y))

def create_emumadz_icon(save_path=None, figsize=(6, 6), dpi=300, 
                        dot_colors=None, dot_sizes=None):
    """
    Create the EMUMADZ icon with customizable colors for each dot
    
    Args:
        dot_colors: List of 9 colors for dots (positions 1-9, left-to-right, top-to-bottom)
        dot_sizes: List of 9 sizes for dots
    """
    
    # Default colors if not provided (positions 1-9: top-left to bottom-right)
    if dot_colors is None:
        dot_colors = ['#0400ff', '#0400ff', '#0400ff',
                      "#757575", '#f39c12', '#757575',
                      "#00ff0d", '#00ff0d', '#00ff0d']
    
    # Default sizes if not provided
    if dot_sizes is None:
        dot_sizes = [400, 400, 400,  # Row 1
                     400, 800, 400,  # Row 2 (center larger)
                     400, 400, 400]  # Row 3
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    
    # Create hexagon boundary
    hex_vertices = create_hexagon(radius=1)
    hexagon = Polygon(hex_vertices, fill=False, edgecolor='#2c3e50', linewidth=3)
    ax.add_patch(hexagon)
    
    # Create 3x3 grid points with custom colors
    x_coords = [-0.4, 0, 0.4]
    y_coords = [0.4, 0, -0.4]  # Top to bottom
    
    # Z pattern connections
    ax.plot([-0.4, 0.4], [0.4, 0.4], color='#34495e', linewidth=2, alpha=0.6)  # Top
    ax.plot([0.4, -0.4], [0.4, -0.4], color='#34495e', linewidth=2, alpha=0.6)  # Diagonal
    ax.plot([-0.4, 0.4], [-0.4, -0.4], color='#34495e', linewidth=2, alpha=0.6)  # Bottom
    
    # Draw grid points with individual colors
    point_index = 0
    for i, y in enumerate(y_coords):
        for j, x in enumerate(x_coords):
            ax.scatter(x, y, s=dot_sizes[point_index], c=dot_colors[point_index], 
                      zorder=5, edgecolors='white', linewidth=1)
            point_index += 1
    
    # Set equal aspect ratio and remove axes
    ax.set_aspect('equal')
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.axis('off')
    
    # Remove margins
    plt.tight_layout()
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
    
    # Save if path provided
    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight', 
                    facecolor='white', edgecolor='none')
        print(f"Icon saved to: {save_path}")
    return fig, ax

if __name__ == "__main__":
    # Generate standard icon
    print("Creating standard EMUMADZ icon...")
    create_emumadz_icon(save_path="emumadz_icon.pdf")
    create_emumadz_icon(save_path="emumadz_icon.png")
