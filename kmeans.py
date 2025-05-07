import math

# Maximum number of iterations for the K-Means algorithm
MAX_ITER = 300

# Threshold for convergence
EPSILON = 0.0001

def EUCLIDEAN_DISTANCE(list1, list2): 
    """
    Calculates the Euclidean distance between two points in multi-dimensional space.

    Parameters:
    list1 (list of float): Coordinates of the first point.
    list2 (list of float): Coordinates of the second point.

    Returns:
    float: Euclidean distance between the two points.
    """
    sum_of_squares = sum((x - y) ** 2 for x, y in zip(list1, list2)) # Calculate the sum of squared differences
    distance = math.sqrt(sum_of_squares) # Compute the Euclidean distance
    return distance 


def FIND_CLOSEST_CENTROID(point , centroids):
  """
    Finds the index of the closest centroid to a given point.

    Parameters:
    point (list of float): The data point.
    centroids (list of list of float): The current centroid coordinates.

    Returns:
    int: Index of the closest centroid.
    """
  
  min_distance = float('inf')
  index_of_closest = 0
  for index_of_centroid in range(len(centroids)):
    distance = EUCLIDEAN_DISTANCE(point , centroids[index_of_centroid] )
    if distance < min_distance: #Update rule
      min_distance = distance
      index_of_closest = index_of_centroid
  return index_of_closest


def CALC_NEW_CORD(points_in_single_centroid , points , cord_i):
  """
    Calculates the average value of a specific coordinate among a list of points.

    Parameters:
    points_in_single_centroid (list of int): Indices of points assigned to the same cluster.
    points (list of list of float): All data points.
    cord_i (int): Index of the coordinate to average.

    Returns:
    float: Averaged coordinate value.
    """
  sum=0
  for index in points_in_single_centroid:
    sum += points[index][cord_i]
  average = sum / len(points_in_single_centroid)
  return average #float(f"{average:.4f}")



def CALC_NEW_CENTROID(points_in_single_centroid , points , centroid):  
  """
    Calculates the new centroid based on the average of all assigned points.

    Parameters:
    points_in_single_centroid (list of int): Indices of points assigned to this centroid.
    points (list of list of float): All data points.
    centroid (list of float): The current centroid (used as fallback if no points assigned).

    Returns:
    list of float: The new centroid.
    """
  if len(points_in_single_centroid) == 0:
    return centroid

  d = len(centroid)
  new_centroid = []
  for cord_i in range (d):
    new_centroid.append(CALC_NEW_CORD(points_in_single_centroid , points , cord_i))
  
  return new_centroid


def UPDATE_POINTS(points , centroids , points_in_centroids):
    """
    Assigns each point to the closest centroid.

    Parameters:
    points (list of list of float): All data points.
    centroids (list of list of float): Current centroids.
    points_in_centroids (list of list of int): Will be filled with indices of assigned points.
    """
    for index_of_point in range(len(points)):
      index_of_closest_centroid = FIND_CLOSEST_CENTROID(points[index_of_point],centroids)
      points_in_centroids[index_of_closest_centroid].append(index_of_point)



def UPDATE_CENTROIDS(centroids , points_in_centroids , points):
  """
    Updates the centroids based on the mean of assigned points and returns the largest centroid movement.

    Parameters:
    centroids (list of list of float): Current centroids.
    points_in_centroids (list of list of int): Points assigned to each centroid.
    points (list of list of float): All data points.

    Returns:
    float: Largest movement (delta) of any centroid.
    """
  largest_delta = 0
  for index_of_centroid in range (len(points_in_centroids)):
    new_centroid = CALC_NEW_CENTROID( points_in_centroids[index_of_centroid] , points , centroids[index_of_centroid])
    temp_delta = EUCLIDEAN_DISTANCE(centroids[index_of_centroid],new_centroid)
    if temp_delta > largest_delta: #update rule for max delta
      largest_delta = temp_delta
    centroids[index_of_centroid] = new_centroid
  return largest_delta


def KmeansAlgorithm(points, k):
  """
    Performs K-Means clustering on a set of points.

    Parameters:
    points (list of list of float): Data points to cluster.
    k (int): Number of clusters.

    Returns:
    list of int: Cluster label for each point (index corresponds to input point index).
  """
  current_iteration=0
  current_max_centroid_delta=float('inf')
  centroids = []
  for i in range (k): #set the initial centroids to a copy of the first k points.  
    centroids.append(list(points[i]))

  while(current_iteration < MAX_ITER and current_max_centroid_delta > EPSILON):
    #reset the points in each centroid.
    points_in_centroids = [[] for i in range(k)] 
      
    #Associate each point with closest cluster and add it to points_in_centroid
    UPDATE_POINTS(points , centroids , points_in_centroids) 
      
    #update each centroid to new location. returns the largest delta
    current_max_centroid_delta = UPDATE_CENTROIDS(centroids , points_in_centroids , points) 
      
    current_iteration+=1
  labels = []
  for j in range(len(points)):
    labels.append(FIND_CLOSEST_CENTROID(points[j], centroids))
  return labels